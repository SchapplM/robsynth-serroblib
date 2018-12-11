% Zeitableitung der geometrischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRPRRR14
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d4,d5,d6,theta3]';
%
% Output:
% JgD [6x6]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-12-10 18:39
% Revision: bb42a8b95257d9bc83910d26e849f5825122f662 (2018-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD = S6RRPRRR14_jacobigD_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)


JaD_transl = S6RRPRRR14_jacobiaD_transl_3_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JgD_rot = S6RRPRRR14_jacobigD_rot_3_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
