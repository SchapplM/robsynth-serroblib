% Analytische Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S4PRRR2
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4]';
%
% Output:
% Ja [6x4]
%   Analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:27
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S4PRRR2_jacobia_3_sym_varpar(qJ, r_i_i_C, ...
  pkin)

Ja_transl = S4PRRR2_jacobia_transl_3_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Ja_rot = S4PRRR2_jacobia_rot_3_sym_varpar(qJ, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
