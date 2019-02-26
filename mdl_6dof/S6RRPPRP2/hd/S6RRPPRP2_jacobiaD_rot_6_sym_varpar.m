% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP2
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:22
% EndTime: 2019-02-26 21:25:22
% DurationCPUTime: 0.74s
% Computational Cost: add. (2081->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->92)
t138 = qJ(2) + pkin(9);
t136 = sin(t138);
t132 = 0.1e1 / t136 ^ 2;
t137 = cos(t138);
t135 = t137 ^ 2;
t188 = t132 * t135;
t204 = t137 * t188;
t142 = sin(qJ(1));
t163 = 0.1e1 + t188;
t203 = t142 * t163;
t144 = cos(qJ(1));
t159 = qJD(1) * t136 + qJD(5);
t175 = qJD(2) * t137;
t202 = t159 * t142 - t144 * t175;
t182 = t142 * t137;
t126 = atan2(-t182, t136);
t125 = cos(t126);
t124 = sin(t126);
t168 = t124 * t182;
t110 = t125 * t136 - t168;
t107 = 0.1e1 / t110;
t143 = cos(qJ(5));
t181 = t142 * t143;
t141 = sin(qJ(5));
t183 = t141 * t144;
t121 = t136 * t183 + t181;
t117 = 0.1e1 / t121;
t131 = 0.1e1 / t136;
t108 = 0.1e1 / t110 ^ 2;
t118 = 0.1e1 / t121 ^ 2;
t139 = t142 ^ 2;
t129 = t139 * t188 + 0.1e1;
t127 = 0.1e1 / t129;
t201 = t127 - 0.1e1;
t177 = qJD(1) * t144;
t166 = t137 * t177;
t174 = qJD(2) * t142;
t101 = ((t136 * t174 - t166) * t131 + t174 * t188) * t127;
t190 = t125 * t137;
t96 = (-t101 * t142 + qJD(2)) * t190 + (-t166 + (-t101 + t174) * t136) * t124;
t200 = t107 * t108 * t96;
t160 = qJD(5) * t136 + qJD(1);
t155 = t160 * t144;
t105 = t141 * t155 + t143 * t202;
t180 = t143 * t144;
t184 = t141 * t142;
t120 = -t136 * t180 + t184;
t116 = t120 ^ 2;
t115 = t116 * t118 + 0.1e1;
t193 = t118 * t120;
t106 = -t141 * t202 + t143 * t155;
t197 = t106 * t117 * t118;
t199 = 0.1e1 / t115 ^ 2 * (t105 * t193 - t116 * t197);
t198 = t101 * t137;
t196 = t108 * t137;
t195 = t108 * t144;
t189 = t131 * t137;
t153 = qJD(2) * (-t131 * t204 - t189);
t186 = t135 * t142;
t157 = t177 * t186;
t194 = (t132 * t157 + t139 * t153) / t129 ^ 2;
t192 = t120 * t141;
t191 = t124 * t136;
t140 = t144 ^ 2;
t187 = t135 * t140;
t185 = t137 * t144;
t112 = t127 * t203;
t179 = -t112 + t142;
t178 = qJD(1) * t142;
t176 = qJD(2) * t136;
t104 = t108 * t187 + 0.1e1;
t173 = 0.2e1 / t104 ^ 2 * (-t187 * t200 + (-t136 * t140 * t175 - t157) * t108);
t172 = 0.2e1 * t200;
t171 = 0.2e1 * t199;
t170 = -0.2e1 * t194;
t169 = t137 * t194;
t167 = t131 * t186;
t164 = t112 * t142 - 0.1e1;
t162 = t137 * t173;
t161 = 0.2e1 * t120 * t197;
t158 = t127 * t167;
t156 = t163 * t144;
t154 = t117 * t143 + t118 * t192;
t152 = t154 * t144;
t123 = -t136 * t184 + t180;
t122 = t136 * t181 + t183;
t113 = 0.1e1 / t115;
t102 = 0.1e1 / t104;
t100 = (t201 * t137 * t124 + t125 * t158) * t144;
t98 = t142 * t191 + t190 + (-t125 * t182 - t191) * t112;
t97 = t170 * t203 + (qJD(1) * t156 + 0.2e1 * t142 * t153) * t127;
t1 = [0.2e1 * t131 * t144 * t169 + (qJD(2) * t156 + t178 * t189) * t127, t97, 0, 0, 0, 0; (t107 * t162 + (t107 * t176 + (qJD(1) * t100 + t96) * t196) * t102) * t142 + (t108 * t162 * t100 + (-((-t101 * t158 - t201 * t176 - 0.2e1 * t169) * t124 + (t167 * t170 - t198 + (t198 + (-0.2e1 * t137 - t204) * t174) * t127) * t125) * t108 * t185 + (t108 * t176 + t137 * t172) * t100 + (-t107 + ((t139 - t140) * t135 * t131 * t127 * t125 + t201 * t168) * t108) * t137 * qJD(1)) * t102) * t144 (t107 * t136 + t98 * t196) * t144 * t173 + ((t107 * t178 + (qJD(2) * t98 + t96) * t195) * t136 + ((-qJD(2) * t107 + t98 * t172) * t144 + (t98 * t178 + (-(-t112 * t177 - t142 * t97) * t125 - (t179 * qJD(2) + t164 * t101) * t124) * t185) * t108 - ((-t97 + t177) * t124 + (t164 * qJD(2) + t179 * t101) * t125) * t136 * t195) * t137) * t102, 0, 0, 0, 0; (-t117 * t122 + t123 * t193) * t171 + (t123 * t161 + (-t123 * t105 - t122 * t106 + t160 * t120 * t181 - (-t137 * t174 - t159 * t144) * t192) * t118 + (t159 * t180 + (-t160 * t141 + t143 * t175) * t142) * t117) * t113, t137 * t152 * t171 + (t152 * t176 + (t154 * t178 + ((qJD(5) * t117 + t161) * t141 + (-t105 * t141 + (-qJD(5) * t120 + t106) * t143) * t118) * t144) * t137) * t113, 0, 0, -0.2e1 * t199 + 0.2e1 * (t105 * t113 * t118 + (-t113 * t197 - t118 * t199) * t120) * t120, 0;];
JaD_rot  = t1;
