% Geometrische Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4RRPP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
%
% Output:
% Jg [3x4]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S4RRPP1_jacobig_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S4RRPP1_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S4RRPP1_jacobig_rot_3_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
