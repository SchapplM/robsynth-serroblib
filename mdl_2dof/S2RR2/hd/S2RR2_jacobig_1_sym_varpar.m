% Geometrische Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S2RR2
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
%
% Output:
% Jg [3x2]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S2RR2_jacobig_1_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S2RR2_jacobia_transl_1_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S2RR2_jacobig_rot_1_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];