//
// E kompleksie vektori
//
// Matricaa U transformaacijai no Dekarta uz ciklisko koordinaatu sisteemu (Auzinsh, Budker, Rochester (3.9), kuru njem kompleksi saistiitu saskanjaa ar (3.15), p. 32), matrica nav atkariiga no E vektora izveeles, taapeec to defineejam globaali
gsl_matrix_complex * U = gsl_matrix_complex_calloc(3, 3);

void init_u(){ // Saliek pareizaas veertiibas U matricaa
    gsl_matrix_complex_set(U, 0, 0, gsl_complex_rect(-1/sqrt(2), 0));
    gsl_matrix_complex_set(U, 0, 1, gsl_complex_rect(0, 1/sqrt(2)));
    gsl_matrix_complex_set(U, 1, 2, gsl_complex_rect(1, 0));
    gsl_matrix_complex_set(U, 2, 0, gsl_complex_rect(1/sqrt(2), 0));
    gsl_matrix_complex_set(U, 2, 1, gsl_complex_rect(0, 1/sqrt(2)));
}

gsl_complex E(int q, int pol_u, double teta_u, double fii_u){
/*
    Shiis funkcijas interfeiss nav mainiits, ir taads pats kaa ieprieksh, tomeer izsaukt to katru reizi kad ievajagaas kaadu E komponenti sanaak neefektiivi, jo katru reizi tiek iziets pils transformaaciju celjsh
*/
    if(abs(q) < 1.5) // ja q par lielu, tad nav jeegas neko dariit
    {
        // inicializeejam saakotneejo vektoru uz nulleem, peec tam saliek veertiibas
        gsl_vector_complex * initial = gsl_vector_complex_calloc(3);
        if(pol_u == 1) // Pa kreisi polarizeets 1/sqrt(2)(x + iy)
        {
            gsl_vector_complex_set(initial, 0, gsl_complex_rect(1/sqrt(2), 0));
            gsl_vector_complex_set(initial, 1, gsl_complex_rect(0, 1/sqrt(2)));
        }
        else if(pol_u == 0) // Lineaars z
        {
            gsl_vector_complex_set(initial, 2, gsl_complex_rect(1, 0));
        }
        else if(pol_u == -1) // Pa labi polarizeets 1/sqrt(2)(x - iy)
        {
            gsl_vector_complex_set(initial, 0, gsl_complex_rect(1/sqrt(2), 0));
            gsl_vector_complex_set(initial, 1, gsl_complex_rect(0, -1/sqrt(2)));
        }
        
        // Pagrieziena matrica ap y asi par teta
        gsl_matrix_complex * R1 = gsl_matrix_complex_calloc(3, 3);
        gsl_matrix_complex_set(R1, 0, 0, gsl_complex_rect(cos(teta_u), 0));
        gsl_matrix_complex_set(R1, 0, 2, gsl_complex_rect(sin(teta_u), 0));
        gsl_matrix_complex_set(R1, 1, 1, gsl_complex_rect(1, 0));
        gsl_matrix_complex_set(R1, 2, 0, gsl_complex_rect(-sin(teta_u), 0));
        gsl_matrix_complex_set(R1, 2, 2, gsl_complex_rect(cos(teta_u), 0));
        
        // Pagrieziena matrica ap z asi par fii
        gsl_matrix_complex * R2 = gsl_matrix_complex_calloc(3, 3);
        gsl_matrix_complex_set(R2, 0, 0, gsl_complex_rect(cos(fii_u), 0));
        gsl_matrix_complex_set(R2, 0, 1, gsl_complex_rect(-sin(fii_u), 0));
        gsl_matrix_complex_set(R2, 1, 0, gsl_complex_rect(sin(fii_u), 0));
        gsl_matrix_complex_set(R2, 1, 1, gsl_complex_rect(cos(fii_u), 0));
        gsl_matrix_complex_set(R2, 2, 2, gsl_complex_rect(1, 0));
        
        // Kopeejaa transformaacijas matrica R = R2 * R1 (transformaacijas matrica iedarbojas uz vektoru no kreisaas puses, taapeec shaada reizinaataaju seciiba)
        gsl_matrix_complex * R = gsl_matrix_complex_calloc(3, 3);
        gsl_blas_zgemm(CblasNoTrans, CblasNoTrans, gsl_complex_rect(1, 0), R2, R1, gsl_complex_rect(0, 0), R);
        
        // Dekarta vektora transformaacija
        gsl_vector_complex * rotated_Cart = gsl_vector_complex_calloc(3);
        gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1, 0), R, initial, gsl_complex_rect(0, 0), rotated_Cart);
        
        // Paareja no Dekarta uz cikliskajaam koordinaataam
        gsl_vector_complex * cyclic = gsl_vector_complex_calloc(3);
        gsl_blas_zgemv(CblasNoTrans, gsl_complex_rect(1, 0), U, rotated_Cart, gsl_complex_rect(0, 0), cyclic);
        
        // Izveelamies veertiibu, ko atgriezt
        gsl_complex retVal = gsl_vector_complex_get(cyclic, -q+1);
        
        // Atbriivojam atminju
        gsl_vector_complex_free(initial);
        gsl_vector_complex_free(rotated_Cart);
        gsl_vector_complex_free(cyclic);

        gsl_matrix_complex_free(R1);
        gsl_matrix_complex_free(R2);
        gsl_matrix_complex_free(R);
        return retVal;
    }
    else
    {
        return gsl_complex_rect(0,0);
    }
}
