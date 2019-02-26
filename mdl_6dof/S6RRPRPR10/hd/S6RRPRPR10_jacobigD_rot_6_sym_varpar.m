% Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JgD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der geometrischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JgD_rot = S6RRPRPR10_jacobigD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR10_jacobigD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR10_jacobigD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR10_jacobigD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobigD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:01
% EndTime: 2019-02-26 21:43:01
% DurationCPUTime: 0.06s
% Computational Cost: add. (33->17), mult. (78->40), div. (0->0), fcn. (80->8), ass. (0->23)
t163 = sin(pkin(6));
t166 = sin(qJ(1));
t181 = t163 * t166;
t168 = cos(qJ(1));
t180 = t163 * t168;
t165 = sin(qJ(2));
t179 = t165 * t166;
t178 = t165 * t168;
t167 = cos(qJ(2));
t177 = t166 * t167;
t176 = t168 * t167;
t175 = qJD(1) * t163;
t162 = pkin(11) + qJ(4);
t161 = cos(t162);
t174 = qJD(2) * t161;
t173 = qJD(2) * t163;
t164 = cos(pkin(6));
t172 = t164 * t176 - t179;
t171 = t164 * t177 + t178;
t170 = t164 * t178 + t177;
t169 = -t164 * t179 + t176;
t160 = sin(t162);
t1 = [0, t168 * t175, 0, t172 * qJD(1) + t169 * qJD(2), 0 (-t169 * t160 + t161 * t181) * qJD(4) - t171 * t174 + (t160 * t180 - t170 * t161) * qJD(1); 0, t166 * t175, 0, t171 * qJD(1) + t170 * qJD(2), 0 (-t170 * t160 - t161 * t180) * qJD(4) + t172 * t174 + (t160 * t181 + t169 * t161) * qJD(1); 0, 0, 0, t165 * t173, 0, t167 * t161 * t173 + (-t160 * t163 * t165 + t161 * t164) * qJD(4);];
JgD_rot  = t1;
