% Analytische Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S3PRP1
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
%
% Output:
% Ja [6x3]
%   Analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:12
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S3PRP1_jacobia_1_sym_varpar(qJ, r_i_i_C, ...
  pkin)

Ja_transl = S3PRP1_jacobia_transl_1_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Ja_rot = S3PRP1_jacobia_rot_1_sym_varpar(qJ, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
