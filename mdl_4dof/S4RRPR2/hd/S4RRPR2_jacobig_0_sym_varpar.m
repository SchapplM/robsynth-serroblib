% Geometrische Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S4RRPR2
% Use Code from Maple symbolic Code Generation
%
% Geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorgeschwindigkeit und Geschw. der verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% Jg [3x4]
%   Geometrische Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S4RRPR2_jacobig_0_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S4RRPR2_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S4RRPR2_jacobig_rot_0_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
