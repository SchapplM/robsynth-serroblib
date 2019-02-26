% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:16
% EndTime: 2019-02-26 20:26:17
% DurationCPUTime: 0.73s
% Computational Cost: add. (3055->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->94)
t142 = pkin(10) + qJ(4);
t138 = sin(t142);
t132 = 0.1e1 / t138 ^ 2;
t140 = cos(t142);
t136 = t140 ^ 2;
t192 = t132 * t136;
t207 = t140 * t192;
t143 = qJ(1) + pkin(9);
t139 = sin(t143);
t165 = 0.1e1 + t192;
t206 = t139 * t165;
t144 = sin(qJ(6));
t145 = cos(qJ(6));
t162 = qJD(6) * t138 + qJD(1);
t177 = qJD(4) * t140;
t205 = t162 * t144 - t145 * t177;
t204 = t144 * t177 + t162 * t145;
t189 = t139 * t140;
t126 = atan2(-t189, t138);
t125 = cos(t126);
t124 = sin(t126);
t171 = t124 * t189;
t112 = t125 * t138 - t171;
t108 = 0.1e1 / t112;
t141 = cos(t143);
t185 = t141 * t144;
t187 = t139 * t145;
t121 = t138 * t185 + t187;
t117 = 0.1e1 / t121;
t131 = 0.1e1 / t138;
t109 = 0.1e1 / t112 ^ 2;
t118 = 0.1e1 / t121 ^ 2;
t134 = t139 ^ 2;
t129 = t134 * t192 + 0.1e1;
t127 = 0.1e1 / t129;
t203 = t127 - 0.1e1;
t180 = qJD(1) * t141;
t169 = t140 * t180;
t178 = qJD(4) * t139;
t101 = ((t138 * t178 - t169) * t131 + t178 * t192) * t127;
t193 = t125 * t140;
t96 = (-t101 * t139 + qJD(4)) * t193 + (-t169 + (-t101 + t178) * t138) * t124;
t202 = t108 * t109 * t96;
t161 = qJD(1) * t138 + qJD(6);
t156 = t161 * t145;
t105 = t139 * t156 + t141 * t205;
t184 = t141 * t145;
t188 = t139 * t144;
t120 = -t138 * t184 + t188;
t116 = t120 ^ 2;
t115 = t116 * t118 + 0.1e1;
t195 = t118 * t120;
t157 = t161 * t144;
t106 = -t139 * t157 + t141 * t204;
t199 = t106 * t117 * t118;
t201 = 0.1e1 / t115 ^ 2 * (t105 * t195 - t116 * t199);
t200 = t101 * t140;
t154 = qJD(4) * (-t140 - t207) * t131;
t190 = t136 * t139;
t159 = t180 * t190;
t198 = (t132 * t159 + t134 * t154) / t129 ^ 2;
t197 = t109 * t140;
t196 = t109 * t141;
t194 = t124 * t138;
t137 = t141 ^ 2;
t191 = t136 * t137;
t186 = t140 * t141;
t111 = t127 * t206;
t183 = -t111 + t139;
t182 = qJD(1) * t139;
t181 = qJD(1) * t140;
t179 = qJD(4) * t138;
t104 = t109 * t191 + 0.1e1;
t176 = 0.2e1 / t104 ^ 2 * (-t191 * t202 + (-t137 * t138 * t177 - t159) * t109);
t175 = 0.2e1 * t202;
t174 = 0.2e1 * t201;
t173 = -0.2e1 * t198;
t172 = t140 * t198;
t170 = t131 * t190;
t166 = t111 * t139 - 0.1e1;
t164 = t140 * t176;
t163 = 0.2e1 * t120 * t199;
t160 = t127 * t170;
t158 = t165 * t141;
t155 = t117 * t145 + t144 * t195;
t153 = t155 * t141;
t123 = -t138 * t188 + t184;
t122 = t138 * t187 + t185;
t113 = 0.1e1 / t115;
t102 = 0.1e1 / t104;
t100 = (t203 * t140 * t124 + t125 * t160) * t141;
t98 = t139 * t194 + t193 + (-t125 * t189 - t194) * t111;
t97 = t173 * t206 + (qJD(1) * t158 + 0.2e1 * t139 * t154) * t127;
t1 = [0.2e1 * t131 * t141 * t172 + (t131 * t139 * t181 + qJD(4) * t158) * t127, 0, 0, t97, 0, 0; (t108 * t164 + (t108 * t179 + (qJD(1) * t100 + t96) * t197) * t102) * t139 + (t109 * t164 * t100 + (-((-t101 * t160 - t203 * t179 - 0.2e1 * t172) * t124 + (t170 * t173 - t200 + (t200 + (-0.2e1 * t140 - t207) * t178) * t127) * t125) * t109 * t186 + (t109 * t179 + t140 * t175) * t100 + (-t108 + ((t134 - t137) * t136 * t131 * t127 * t125 + t203 * t171) * t109) * t181) * t102) * t141, 0, 0 (t108 * t138 + t98 * t197) * t141 * t176 + ((t108 * t182 + (qJD(4) * t98 + t96) * t196) * t138 + ((-qJD(4) * t108 + t98 * t175) * t141 + (t98 * t182 + (-(-t111 * t180 - t139 * t97) * t125 - (t183 * qJD(4) + t166 * t101) * t124) * t186) * t109 - ((-t97 + t180) * t124 + (t166 * qJD(4) + t183 * t101) * t125) * t138 * t196) * t140) * t102, 0, 0; (-t117 * t122 + t123 * t195) * t174 + (t123 * t163 + (-t123 * t105 - t122 * t106 + (t139 * t204 + t141 * t157) * t120) * t118 + (-t205 * t139 + t141 * t156) * t117) * t113, 0, 0, t140 * t153 * t174 + (t153 * t179 + (t155 * t182 + ((qJD(6) * t117 + t163) * t144 + (-t105 * t144 + (-qJD(6) * t120 + t106) * t145) * t118) * t141) * t140) * t113, 0, -0.2e1 * t201 + 0.2e1 * (t105 * t113 * t118 + (-t113 * t199 - t118 * t201) * t120) * t120;];
JaD_rot  = t1;
