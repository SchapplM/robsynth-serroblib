% Zeitableitung der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPP3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1]';
%
% Output:
% JgD [6x6]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:22
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function JgD = S6PRRRPP3_jacobigD_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)


JaD_transl = S6PRRRPP3_jacobiaD_transl_6_floatb_twist_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JgD_rot = S6PRRRPP3_jacobigD_rot_6_floatb_twist_sym_varpar(qJ, qJD, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
