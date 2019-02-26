% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPR3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:47
% EndTime: 2019-02-26 22:31:48
% DurationCPUTime: 0.60s
% Computational Cost: add. (6718->70), mult. (3877->159), div. (936->14), fcn. (4495->7), ass. (0->76)
t129 = sin(qJ(1));
t123 = t129 ^ 2;
t122 = qJ(2) + qJ(3) + qJ(4);
t119 = sin(t122);
t114 = t119 ^ 2;
t120 = cos(t122);
t117 = 0.1e1 / t120 ^ 2;
t167 = t114 * t117;
t109 = t123 * t167 + 0.1e1;
t113 = t119 * t114;
t115 = t120 ^ 2;
t116 = 0.1e1 / t120;
t121 = qJD(2) + qJD(3) + qJD(4);
t165 = t116 * t119;
t137 = t121 * (t113 * t116 / t115 + t165);
t130 = cos(qJ(1));
t157 = qJD(1) * t130;
t150 = t129 * t157;
t171 = 0.1e1 / t109 ^ 2 * (t123 * t137 + t150 * t167);
t178 = -0.2e1 * t171;
t162 = t121 * t129;
t107 = 0.1e1 / t109;
t149 = t119 * t157;
t151 = t117 * t162;
t93 = (-(-t120 * t162 - t149) * t116 + t114 * t151) * t107;
t143 = t93 - t162;
t161 = 0.1e1 / t129 * t130;
t147 = 0.1e1 + t167;
t177 = t129 * t147;
t128 = t130 ^ 2;
t176 = qJD(1) * (0.1e1 / t123 * t128 + 0.1e1) * t161;
t159 = t129 * t119;
t106 = atan2(-t159, -t120);
t105 = cos(t106);
t104 = sin(t106);
t152 = t104 * t159;
t101 = -t105 * t120 - t152;
t98 = 0.1e1 / t101;
t99 = 0.1e1 / t101 ^ 2;
t175 = t107 - 0.1e1;
t163 = t120 * t121;
t141 = t119 * t128 * t163;
t144 = -t129 * t93 + t121;
t168 = t105 * t119;
t88 = t144 * t168 + (t143 * t120 - t149) * t104;
t172 = t98 * t99 * t88;
t96 = t128 * t114 * t99 + 0.1e1;
t174 = (t99 * t141 + (-t128 * t172 - t99 * t150) * t114) / t96 ^ 2;
t94 = 0.1e1 / t96;
t173 = t94 * t99;
t170 = t130 * t99;
t169 = t104 * t129;
t166 = t114 * t129;
t164 = t119 * t130;
t125 = 0.1e1 / t129 ^ 2;
t160 = t125 * t128;
t158 = qJD(1) * t129;
t156 = 0.2e1 * t172;
t155 = t98 * t174;
t112 = t115 * t160 + 0.1e1;
t154 = 0.2e1 * (-t115 * t176 - t125 * t141) / t112 ^ 2;
t153 = t94 * t163;
t148 = 0.2e1 * t99 * t174;
t146 = 0.1e1 + t160;
t145 = t116 * t178;
t142 = t105 * t107 * t114 * t116;
t140 = t147 * t130;
t139 = t146 * t119;
t110 = 0.1e1 / t112;
t102 = t107 * t177;
t92 = (t175 * t119 * t104 - t129 * t142) * t130;
t91 = t119 * t154 * t161 + (qJD(1) * t139 - t161 * t163) * t110;
t90 = -t120 * t169 + t168 + (t104 * t120 - t105 * t159) * t102;
t89 = t177 * t178 + (qJD(1) * t140 + 0.2e1 * t129 * t137) * t107;
t86 = (-t98 * t94 * t158 + (-0.2e1 * t155 + (-t121 * t90 - t88) * t173) * t130) * t120 + (t90 * t130 * t148 + (-t130 * t121 * t98 - (-t105 * t129 * t89 - t143 * t104 + (-t104 * t121 - t105 * t157 + t169 * t93) * t102) * t99 * t164 + (t130 * t156 + t99 * t158) * t90 - ((t89 - t157) * t104 + (t143 * t102 + t144) * t105) * t120 * t170) * t94) * t119;
t1 = [t145 * t164 + (t121 * t140 - t158 * t165) * t107, t89, t89, t89, 0, 0; (-t98 * t153 + (0.2e1 * t155 + (qJD(1) * t92 + t88) * t173) * t119) * t129 + (-t92 * t99 * t153 + (t92 * t148 + (t92 * t156 + ((0.2e1 * t119 * t171 + t163 + (-t116 * t93 * t166 - t163) * t107) * t104 + (t145 * t166 + t119 * t93 + (t113 * t151 - (t93 - 0.2e1 * t162) * t119) * t107) * t105) * t170) * t94) * t119 + (-t98 + (-(t123 - t128) * t142 + t175 * t152) * t99) * t119 * t94 * qJD(1)) * t130, t86, t86, t86, 0, 0; t146 * t120 * t154 + (0.2e1 * t120 * t176 + t121 * t139) * t110, t91, t91, t91, 0, 0;];
JaD_rot  = t1;
