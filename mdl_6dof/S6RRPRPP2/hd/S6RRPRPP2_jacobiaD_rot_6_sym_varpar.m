% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPP2_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:21
% EndTime: 2019-02-26 21:35:22
% DurationCPUTime: 0.73s
% Computational Cost: add. (1907->96), mult. (2545->208), div. (484->12), fcn. (2996->9), ass. (0->94)
t139 = qJ(2) + pkin(9);
t137 = sin(t139);
t133 = t137 ^ 2;
t138 = cos(t139);
t135 = 0.1e1 / t138 ^ 2;
t187 = t133 * t135;
t143 = sin(qJ(1));
t140 = t143 ^ 2;
t129 = t140 * t187 + 0.1e1;
t134 = 0.1e1 / t138;
t184 = t134 * t137;
t203 = t137 * t187;
t153 = qJD(2) * (t134 * t203 + t184);
t145 = cos(qJ(1));
t175 = qJD(1) * t145;
t185 = t133 * t143;
t156 = t175 * t185;
t194 = (t135 * t156 + t140 * t153) / t129 ^ 2;
t205 = -0.2e1 * t194;
t142 = sin(qJ(4));
t144 = cos(qJ(4));
t172 = qJD(2) * t145;
t164 = t137 * t172;
t176 = qJD(1) * t143;
t179 = t145 * t142;
t180 = t143 * t144;
t106 = t144 * t164 - qJD(4) * t180 - t142 * t175 + (qJD(4) * t179 + t144 * t176) * t138;
t122 = -t138 * t179 + t180;
t116 = t122 ^ 2;
t178 = t145 * t144;
t181 = t143 * t142;
t123 = t138 * t178 + t181;
t118 = 0.1e1 / t123 ^ 2;
t193 = t116 * t118;
t115 = 0.1e1 + t193;
t158 = qJD(1) * t138 - qJD(4);
t159 = qJD(4) * t138 - qJD(1);
t105 = -t159 * t178 + (t158 * t143 + t164) * t142;
t190 = t118 * t122;
t169 = t105 * t190;
t117 = 0.1e1 / t123;
t119 = t117 * t118;
t192 = t116 * t119;
t204 = 0.1e1 / t115 ^ 2 * (t106 * t192 + t169);
t162 = 0.1e1 + t187;
t202 = t143 * t162;
t182 = t143 * t137;
t128 = atan2(t182, t138);
t125 = cos(t128);
t124 = sin(t128);
t167 = t124 * t182;
t110 = t125 * t138 + t167;
t107 = 0.1e1 / t110;
t108 = 0.1e1 / t110 ^ 2;
t201 = 0.2e1 * t137;
t126 = 0.1e1 / t129;
t200 = t126 - 0.1e1;
t141 = t145 ^ 2;
t186 = t133 * t141;
t104 = t108 * t186 + 0.1e1;
t174 = qJD(2) * t138;
t165 = t137 * t175;
t173 = qJD(2) * t143;
t101 = ((t138 * t173 + t165) * t134 + t173 * t187) * t126;
t188 = t125 * t137;
t96 = (t101 * t143 - qJD(2)) * t188 + (t165 + (-t101 + t173) * t138) * t124;
t198 = t107 * t108 * t96;
t199 = 0.1e1 / t104 ^ 2 * (-t186 * t198 + (t137 * t141 * t174 - t156) * t108);
t197 = t101 * t137;
t196 = t108 * t137;
t195 = t108 * t145;
t191 = t117 * t142;
t189 = t124 * t138;
t183 = t137 * t145;
t112 = t126 * t202;
t177 = t112 - t143;
t171 = -0.2e1 * t198;
t170 = 0.2e1 * t204;
t168 = t122 * t119 * t106;
t166 = t126 * t133 * t134;
t163 = t112 * t143 - 0.1e1;
t161 = -0.2e1 * t137 * t199;
t160 = t134 * t205;
t157 = t143 * t166;
t155 = t162 * t145;
t154 = t144 * t190 + t191;
t121 = -t138 * t180 + t179;
t120 = t138 * t181 + t178;
t113 = 0.1e1 / t115;
t102 = 0.1e1 / t104;
t100 = (-t200 * t137 * t124 + t125 * t157) * t145;
t98 = t143 * t189 - t188 + (t125 * t182 - t189) * t112;
t97 = t202 * t205 + (qJD(1) * t155 + 0.2e1 * t143 * t153) * t126;
t1 = [t160 * t183 + (qJD(2) * t155 - t176 * t184) * t126, t97, 0, 0, 0, 0; (t107 * t161 + (t107 * t174 + (-qJD(1) * t100 - t96) * t196) * t102) * t143 + (t108 * t161 * t100 + (((-t101 * t157 - t200 * t174 + t194 * t201) * t124 + (t160 * t185 + t197 + (-t197 + (t201 + t203) * t173) * t126) * t125) * t108 * t183 + (t108 * t174 + t137 * t171) * t100 + (t107 + ((-t140 + t141) * t125 * t166 + t200 * t167) * t108) * t137 * qJD(1)) * t102) * t145, 0.2e1 * (t107 * t138 - t98 * t196) * t145 * t199 + ((t107 * t176 + (qJD(2) * t98 + t96) * t195) * t138 + ((qJD(2) * t107 + t98 * t171) * t145 + (-t98 * t176 + ((t112 * t175 + t143 * t97) * t125 + (t177 * qJD(2) - t163 * t101) * t124) * t183) * t108 + ((-t97 + t175) * t124 + (t163 * qJD(2) - t177 * t101) * t125) * t138 * t195) * t137) * t102, 0, 0, 0, 0; (-t117 * t120 + t121 * t190) * t170 + (-0.2e1 * t121 * t168 + t159 * t117 * t180 + (-t137 * t173 + t158 * t145) * t191 + (-t121 * t105 + t120 * t106 + (t158 * t178 - (qJD(2) * t137 * t144 + t159 * t142) * t143) * t122) * t118) * t113, -0.2e1 * t154 * t183 * t204 + (t154 * t138 * t172 + (-t154 * t176 + ((qJD(4) * t117 + 0.2e1 * t168) * t144 + (t105 * t144 + (-qJD(4) * t122 + t106) * t142) * t118) * t145) * t137) * t113, 0 (t117 * t123 + t193) * t170 + (-0.2e1 * t169 + (-t118 * t123 + t117 - 0.2e1 * t192) * t106) * t113, 0, 0;];
JaD_rot  = t1;
