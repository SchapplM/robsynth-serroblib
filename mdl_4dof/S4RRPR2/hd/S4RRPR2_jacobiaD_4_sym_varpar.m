% Zeitableitung der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S4RRPR2
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
% qJD [4x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
%
% Output:
% JaD [6x4]
%   Zeitableitung der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 18:16
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD = S4RRPR2_jacobiaD_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)

JaD_transl = S4RRPR2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JaD_rot = S4RRPR2_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin);

JaD = [JaD_transl; JaD_rot];
