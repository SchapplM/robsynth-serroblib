% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6PRPRRR6_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:55
% EndTime: 2019-02-26 19:56:55
% DurationCPUTime: 0.04s
% Computational Cost: add. (22->11), mult. (84->29), div. (0->0), fcn. (88->8), ass. (0->19)
t155 = sin(pkin(6));
t160 = cos(qJ(4));
t168 = t155 * t160;
t161 = cos(qJ(2));
t167 = t155 * t161;
t157 = cos(pkin(6));
t159 = sin(qJ(2));
t166 = t157 * t159;
t165 = t157 * t161;
t164 = qJD(2) * t160;
t154 = sin(pkin(11));
t156 = cos(pkin(11));
t163 = -t154 * t159 + t156 * t165;
t162 = t154 * t165 + t156 * t159;
t158 = sin(qJ(4));
t153 = -t155 * t159 * t164 + (t157 * t160 - t158 * t167) * qJD(4);
t152 = (-t156 * t168 - t163 * t158) * qJD(4) - (t154 * t161 + t156 * t166) * t164;
t151 = (t154 * t168 + t162 * t158) * qJD(4) - (-t154 * t166 + t156 * t161) * t164;
t1 = [0, 0, 0, -t162 * qJD(2), t151, t151; 0, 0, 0, t163 * qJD(2), t152, t152; 0, 0, 0, qJD(2) * t167, t153, t153;];
JgD_rot  = t1;
