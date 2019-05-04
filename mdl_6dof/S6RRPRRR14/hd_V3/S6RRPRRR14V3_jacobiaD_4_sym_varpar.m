% Zeitableitung der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD [6x6]
%   Zeitableitung der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD = S6RRPRRR14V3_jacobiaD_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)

JaD_transl = S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JaD_rot = S6RRPRRR14V3_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin);

JaD = [JaD_transl; JaD_rot];
