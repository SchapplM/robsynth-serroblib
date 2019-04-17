% Analytische Jacobi-Matrix für Segment Nr. 0 (0=Basis) von
% S3PPP1
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
%   pkin=[a2,a3,theta1]';
%
% Output:
% Ja [6x3]
%   Analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-17 09:48
% Revision: 3acd05283b8979b361f80d69cfa1c98d98241298 (2019-04-12)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function Ja = S3PPP1_jacobia_0_sym_varpar(qJ, r_i_i_C, ...
  pkin)

Ja_transl = S3PPP1_jacobia_transl_0_sym_varpar(qJ, r_i_i_C, ...
  pkin);
Ja_rot = S3PPP1_jacobia_rot_0_sym_varpar(qJ, ...
  pkin);

Ja = [Ja_transl; Ja_rot];
