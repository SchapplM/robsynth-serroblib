% Geometrische Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S3RPP1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1]';
%
% Output:
% Jg [3x3]
%   Geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:14
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Jg = S3RPP1_jacobig_0_sym_varpar(qJ, r_i_i_C, ...
  pkin)


Ja_transl = S3RPP1_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Jg_rot = S3RPP1_jacobig_rot_0_sym_varpar(qJ, ...
  pkin);

Jg = [Ja_transl; Jg_rot];
