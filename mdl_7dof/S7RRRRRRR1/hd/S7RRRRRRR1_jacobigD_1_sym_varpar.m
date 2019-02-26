% Zeitableitung der geometrischen Jacobi-Matrix für Segment Nr. 1 (0=Basis) von
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
%
% Output:
% JgD [6x7]
%   Zeitableitung der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:54
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD = S7RRRRRRR1_jacobigD_1_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)


JaD_transl = S7RRRRRRR1_jacobiaD_transl_1_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin);
JgD_rot = S7RRRRRRR1_jacobigD_rot_1_sym_varpar(qJ, qJD, ...
  pkin);

JgD = [JaD_transl; JgD_rot];
