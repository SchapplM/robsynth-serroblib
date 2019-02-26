% Zeitableitung der geometrischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JgD [6x6]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD = S6RRPRPR7_jacobigD_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)


JaD_transl = S6RRPRPR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JgD_rot = S6RRPRPR7_jacobigD_rot_4_sym_varpar(qJ, qJD, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
