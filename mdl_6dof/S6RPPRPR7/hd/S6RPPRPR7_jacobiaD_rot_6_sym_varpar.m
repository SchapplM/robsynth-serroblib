% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:09
% EndTime: 2019-02-26 20:29:10
% DurationCPUTime: 0.71s
% Computational Cost: add. (2429->92), mult. (2519->204), div. (480->12), fcn. (2968->9), ass. (0->95)
t141 = pkin(9) + qJ(4);
t137 = sin(t141);
t132 = 0.1e1 / t137 ^ 2;
t139 = cos(t141);
t135 = t139 ^ 2;
t189 = t132 * t135;
t145 = cos(qJ(1));
t206 = 0.2e1 * t145;
t205 = t139 * t189;
t182 = t145 * t139;
t126 = atan2(-t182, t137);
t124 = sin(t126);
t125 = cos(t126);
t110 = -t124 * t182 + t125 * t137;
t107 = 0.1e1 / t110;
t140 = pkin(10) + qJ(6);
t136 = sin(t140);
t138 = cos(t140);
t144 = sin(qJ(1));
t184 = t144 * t138;
t121 = t136 * t145 + t137 * t184;
t117 = 0.1e1 / t121;
t131 = 0.1e1 / t137;
t108 = 0.1e1 / t110 ^ 2;
t118 = 0.1e1 / t121 ^ 2;
t143 = t145 ^ 2;
t129 = t143 * t189 + 0.1e1;
t127 = 0.1e1 / t129;
t204 = t127 - 0.1e1;
t142 = t144 ^ 2;
t188 = t135 * t142;
t106 = t108 * t188 + 0.1e1;
t179 = qJD(1) * t145;
t161 = t135 * t144 * t179;
t177 = qJD(4) * t139;
t180 = qJD(1) * t144;
t170 = t139 * t180;
t176 = qJD(4) * t145;
t101 = ((t137 * t176 + t170) * t131 + t176 * t189) * t127;
t192 = t125 * t139;
t96 = (-t101 * t145 + qJD(4)) * t192 + (t170 + (-t101 + t176) * t137) * t124;
t202 = t107 * t108 * t96;
t203 = 0.1e1 / t106 ^ 2 * (-t188 * t202 + (-t137 * t142 * t177 + t161) * t108);
t164 = qJD(1) * t137 + qJD(6);
t153 = t144 * t177 + t164 * t145;
t165 = qJD(6) * t137 + qJD(1);
t158 = t138 * t165;
t102 = t153 * t136 + t144 * t158;
t183 = t145 * t138;
t185 = t144 * t136;
t120 = t137 * t185 - t183;
t116 = t120 ^ 2;
t115 = t116 * t118 + 0.1e1;
t194 = t118 * t120;
t159 = t136 * t165;
t103 = t153 * t138 - t144 * t159;
t199 = t103 * t117 * t118;
t201 = 0.1e1 / t115 ^ 2 * (t102 * t194 - t116 * t199);
t200 = t101 * t139;
t198 = t108 * t139;
t197 = t108 * t144;
t190 = t131 * t139;
t156 = qJD(4) * (-t131 * t205 - t190);
t196 = (-t132 * t161 + t143 * t156) / t129 ^ 2;
t195 = t117 * t136;
t193 = t120 * t138;
t191 = t131 * t135;
t187 = t137 * t145;
t186 = t139 * t144;
t168 = 0.1e1 + t189;
t114 = t168 * t145 * t127;
t181 = -t114 + t145;
t178 = qJD(4) * t137;
t175 = -0.2e1 * t202;
t174 = 0.2e1 * t201;
t173 = t139 * t203;
t172 = t139 * t196;
t171 = t127 * t191;
t169 = t114 * t145 - 0.1e1;
t167 = 0.2e1 * t120 * t199;
t166 = t196 * t206;
t163 = t145 * t171;
t162 = t204 * t139 * t124;
t160 = t168 * t144;
t157 = t118 * t193 - t195;
t155 = t157 * t144;
t154 = t139 * t176 - t164 * t144;
t123 = t137 * t183 - t185;
t122 = t136 * t187 + t184;
t112 = 0.1e1 / t115;
t104 = 0.1e1 / t106;
t100 = (-t125 * t163 - t162) * t144;
t99 = t124 * t187 + t192 + (-t124 * t137 - t125 * t182) * t114;
t97 = -t168 * t166 + (-qJD(1) * t160 + t156 * t206) * t127;
t1 = [-0.2e1 * t144 * t131 * t172 + (-qJD(4) * t160 + t179 * t190) * t127, 0, 0, t97, 0, 0; (0.2e1 * t107 * t173 + (t107 * t178 + (qJD(1) * t100 + t96) * t198) * t104) * t145 + (-0.2e1 * t108 * t173 * t100 + (((t101 * t163 + t204 * t178 + 0.2e1 * t172) * t124 + (t166 * t191 + t200 + (-t200 + (0.2e1 * t139 + t205) * t176) * t127) * t125) * t108 * t186 + (-t108 * t178 + t139 * t175) * t100 + (t107 + ((t142 - t143) * t125 * t171 - t145 * t162) * t108) * t139 * qJD(1)) * t104) * t144, 0, 0, 0.2e1 * (-t107 * t137 - t99 * t198) * t144 * t203 + ((t107 * t179 + (-qJD(4) * t99 - t96) * t197) * t137 + ((qJD(4) * t107 + t99 * t175) * t144 + (t99 * t179 + ((t114 * t180 - t145 * t97) * t125 + (t181 * qJD(4) + t169 * t101) * t124) * t186) * t108 + ((-t97 - t180) * t124 + (t169 * qJD(4) + t181 * t101) * t125) * t137 * t197) * t139) * t104, 0, 0; (-t117 * t122 + t123 * t194) * t174 + (t123 * t167 + t145 * t117 * t158 + t154 * t195 + (t145 * t120 * t159 - t123 * t102 - t122 * t103 - t154 * t193) * t118) * t112, 0, 0, t139 * t155 * t174 + (t155 * t178 + (-t157 * t179 + ((qJD(6) * t117 + t167) * t138 + (-t102 * t138 + (qJD(6) * t120 - t103) * t136) * t118) * t144) * t139) * t112, 0, -0.2e1 * t201 + 0.2e1 * (t102 * t112 * t118 + (-t112 * t199 - t118 * t201) * t120) * t120;];
JaD_rot  = t1;
