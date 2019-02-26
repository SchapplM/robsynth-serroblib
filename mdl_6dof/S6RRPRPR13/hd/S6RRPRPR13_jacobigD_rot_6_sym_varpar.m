% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR13_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:44:50
% EndTime: 2019-02-26 21:44:50
% DurationCPUTime: 0.06s
% Computational Cost: add. (22->16), mult. (78->39), div. (0->0), fcn. (80->8), ass. (0->22)
t159 = sin(pkin(6));
t163 = sin(qJ(1));
t179 = t159 * t163;
t165 = cos(qJ(2));
t178 = t159 * t165;
t166 = cos(qJ(1));
t177 = t159 * t166;
t162 = sin(qJ(2));
t176 = t163 * t162;
t175 = t163 * t165;
t174 = t165 * t166;
t173 = t166 * t162;
t172 = qJD(1) * t159;
t164 = cos(qJ(4));
t171 = qJD(2) * t164;
t160 = cos(pkin(6));
t170 = t160 * t174 - t176;
t169 = t160 * t175 + t173;
t168 = t160 * t173 + t175;
t167 = -t160 * t176 + t174;
t161 = sin(qJ(4));
t1 = [0, t166 * t172, 0, -t168 * qJD(1) - t169 * qJD(2), 0 (t169 * t161 + t164 * t179) * qJD(4) - t167 * t171 + (t161 * t177 - t170 * t164) * qJD(1); 0, t163 * t172, 0, t167 * qJD(1) + t170 * qJD(2), 0 (-t170 * t161 - t164 * t177) * qJD(4) - t168 * t171 + (t161 * t179 - t169 * t164) * qJD(1); 0, 0, 0, qJD(2) * t178, 0, -t159 * t162 * t171 + (t160 * t164 - t161 * t178) * qJD(4);];
JgD_rot  = t1;
