% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPRR7
% Use Code from Maple symbolic Code Generation
%
% geometrische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRRPRR7_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR7_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR7_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRR7_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:07:37
% EndTime: 2019-02-26 20:07:37
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->10), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t156 = sin(pkin(6));
t160 = sin(qJ(2));
t169 = t156 * t160;
t161 = cos(qJ(3));
t168 = t156 * t161;
t158 = cos(pkin(6));
t167 = t158 * t160;
t162 = cos(qJ(2));
t166 = t158 * t162;
t165 = qJD(2) * t161;
t155 = sin(pkin(11));
t157 = cos(pkin(11));
t164 = t155 * t162 + t157 * t167;
t163 = -t155 * t167 + t157 * t162;
t159 = sin(qJ(3));
t154 = t156 * t162 * t165 + (t158 * t161 - t159 * t169) * qJD(3);
t153 = (t155 * t168 - t163 * t159) * qJD(3) + (-t155 * t166 - t157 * t160) * t165;
t152 = (-t157 * t168 - t164 * t159) * qJD(3) + (-t155 * t160 + t157 * t166) * t165;
t1 = [0, 0, t163 * qJD(2), 0, t153, t153; 0, 0, t164 * qJD(2), 0, t152, t152; 0, 0, qJD(2) * t169, 0, t154, t154;];
JgD_rot  = t1;
