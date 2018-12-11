% Geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% Jg [3x6]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Jg = S6RRPRRR14_jacobig_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S6RRPRRR14_jacobia_transl_3_floatb_twist_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S6RRPRRR14_jacobig_rot_3_floatb_twist_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
