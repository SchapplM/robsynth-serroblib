% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
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

function JaD_rot = S6RRPPRP2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:22
% EndTime: 2019-02-26 21:25:22
% DurationCPUTime: 0.72s
% Computational Cost: add. (2081->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->90)
t138 = sin(qJ(1));
t140 = cos(qJ(1));
t134 = qJ(2) + pkin(9);
t132 = sin(t134);
t155 = qJD(1) * t132 + qJD(5);
t133 = cos(t134);
t175 = qJD(2) * t133;
t200 = t155 * t138 - t140 * t175;
t181 = t138 * t133;
t122 = atan2(-t181, t132);
t121 = cos(t122);
t120 = sin(t122);
t168 = t120 * t181;
t106 = t121 * t132 - t168;
t103 = 0.1e1 / t106;
t139 = cos(qJ(5));
t180 = t138 * t139;
t137 = sin(qJ(5));
t182 = t137 * t140;
t117 = t132 * t182 + t180;
t113 = 0.1e1 / t117;
t127 = 0.1e1 / t132;
t104 = 0.1e1 / t106 ^ 2;
t114 = 0.1e1 / t117 ^ 2;
t128 = 0.1e1 / t132 ^ 2;
t135 = t138 ^ 2;
t131 = t133 ^ 2;
t186 = t128 * t131;
t125 = t135 * t186 + 0.1e1;
t123 = 0.1e1 / t125;
t199 = t123 - 0.1e1;
t136 = t140 ^ 2;
t185 = t131 * t136;
t100 = t104 * t185 + 0.1e1;
t98 = 0.1e1 / t100;
t198 = t104 * t98;
t174 = qJD(2) * t138;
t165 = t128 * t174;
t177 = qJD(1) * t140;
t166 = t133 * t177;
t97 = ((t132 * t174 - t166) * t127 + t131 * t165) * t123;
t157 = -t97 + t174;
t158 = -t138 * t97 + qJD(2);
t189 = t121 * t133;
t92 = t158 * t189 + (t157 * t132 - t166) * t120;
t197 = t103 * t104 * t92;
t156 = qJD(5) * t132 + qJD(1);
t151 = t156 * t140;
t101 = t137 * t151 + t139 * t200;
t179 = t139 * t140;
t183 = t137 * t138;
t116 = -t132 * t179 + t183;
t112 = t116 ^ 2;
t111 = t112 * t114 + 0.1e1;
t192 = t114 * t116;
t102 = -t137 * t200 + t139 * t151;
t194 = t102 * t113 * t114;
t196 = 0.1e1 / t111 ^ 2 * (t101 * t192 - t112 * t194);
t191 = t116 * t137;
t190 = t120 * t132;
t188 = t127 * t131;
t187 = t127 * t133;
t178 = qJD(1) * t138;
t176 = qJD(2) * t132;
t153 = t131 * t138 * t177;
t173 = 0.2e1 * (-t185 * t197 + (-t132 * t136 * t175 - t153) * t104) / t100 ^ 2;
t172 = 0.2e1 * t197;
t171 = 0.2e1 * t196;
t130 = t133 * t131;
t149 = qJD(2) * (-t127 * t128 * t130 - t187);
t170 = 0.2e1 * (t128 * t153 + t135 * t149) / t125 ^ 2;
t169 = t98 * t176;
t167 = t123 * t188;
t163 = 0.1e1 + t186;
t162 = t103 * t173;
t161 = 0.2e1 * t116 * t194;
t160 = t133 * t170;
t159 = t138 * t170;
t154 = t138 * t167;
t152 = t163 * t140;
t150 = t113 * t139 + t114 * t191;
t148 = t150 * t140;
t119 = -t132 * t183 + t179;
t118 = t132 * t180 + t182;
t109 = 0.1e1 / t111;
t108 = t163 * t138 * t123;
t96 = (t199 * t133 * t120 + t121 * t154) * t140;
t94 = t138 * t190 + t189 + (-t121 * t181 - t190) * t108;
t93 = -t163 * t159 + (qJD(1) * t152 + 0.2e1 * t138 * t149) * t123;
t1 = [t127 * t140 * t160 + (qJD(2) * t152 + t178 * t187) * t123, t93, 0, 0, 0, 0; (t103 * t169 + (t162 + (qJD(1) * t96 + t92) * t198) * t133) * t138 + ((t96 * t169 + (t96 * t173 + ((t97 * t154 + t199 * t176 + t160) * t120 + (t159 * t188 + t133 * t97 + (t130 * t165 - (t97 - 0.2e1 * t174) * t133) * t123) * t121) * t98 * t140) * t133) * t104 + (t96 * t172 + (-t103 + ((t135 - t136) * t121 * t167 + t199 * t168) * t104) * qJD(1)) * t133 * t98) * t140 (t103 * t98 * t178 + (t162 + (qJD(2) * t94 + t92) * t198) * t140) * t132 + (((-qJD(2) * t103 + t94 * t172) * t140 + (t94 * t178 + (-(-t108 * t177 - t138 * t93) * t121 - ((t108 * t138 - 0.1e1) * t97 + (-t108 + t138) * qJD(2)) * t120) * t133 * t140) * t104) * t98 + (t94 * t173 - ((-t93 + t177) * t120 + (t157 * t108 - t158) * t121) * t98 * t132) * t104 * t140) * t133, 0, 0, 0, 0; (-t113 * t118 + t119 * t192) * t171 + (t119 * t161 + (-t119 * t101 - t118 * t102 + t156 * t116 * t180 - (-t133 * t174 - t155 * t140) * t191) * t114 + (t155 * t179 + (-t156 * t137 + t139 * t175) * t138) * t113) * t109, t133 * t148 * t171 + (t148 * t176 + (t150 * t178 + ((qJD(5) * t113 + t161) * t137 + (-t101 * t137 + (-qJD(5) * t116 + t102) * t139) * t114) * t140) * t133) * t109, 0, 0, -0.2e1 * t196 + 0.2e1 * (t101 * t109 * t114 + (-t109 * t194 - t114 * t196) * t116) * t116, 0;];
JaD_rot  = t1;
