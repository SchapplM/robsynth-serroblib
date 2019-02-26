% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:38
% EndTime: 2019-02-26 20:58:39
% DurationCPUTime: 0.68s
% Computational Cost: add. (1907->96), mult. (2545->213), div. (484->12), fcn. (2996->9), ass. (0->94)
t140 = sin(qJ(1));
t137 = t140 ^ 2;
t136 = pkin(9) + qJ(3);
t134 = sin(t136);
t130 = t134 ^ 2;
t135 = cos(t136);
t132 = 0.1e1 / t135 ^ 2;
t187 = t130 * t132;
t126 = t137 * t187 + 0.1e1;
t129 = t134 * t130;
t131 = 0.1e1 / t135;
t184 = t131 * t134;
t150 = qJD(3) * (t129 * t131 * t132 + t184);
t142 = cos(qJ(1));
t176 = qJD(1) * t142;
t185 = t130 * t140;
t153 = t176 * t185;
t194 = (t132 * t153 + t137 * t150) / t126 ^ 2;
t203 = -0.2e1 * t194;
t139 = sin(qJ(4));
t141 = cos(qJ(4));
t173 = qJD(3) * t142;
t161 = t134 * t173;
t177 = qJD(1) * t140;
t179 = t140 * t141;
t182 = t139 * t142;
t103 = t141 * t161 - qJD(4) * t179 - t139 * t176 + (qJD(4) * t182 + t141 * t177) * t135;
t119 = -t135 * t182 + t179;
t113 = t119 ^ 2;
t178 = t141 * t142;
t180 = t140 * t139;
t120 = t135 * t178 + t180;
t115 = 0.1e1 / t120 ^ 2;
t193 = t113 * t115;
t112 = 0.1e1 + t193;
t155 = qJD(1) * t135 - qJD(4);
t156 = qJD(4) * t135 - qJD(1);
t102 = -t156 * t178 + (t155 * t140 + t161) * t139;
t190 = t115 * t119;
t167 = t102 * t190;
t114 = 0.1e1 / t120;
t116 = t114 * t115;
t192 = t113 * t116;
t202 = 0.1e1 / t112 ^ 2 * (t103 * t192 + t167);
t160 = 0.1e1 + t187;
t201 = t140 * t160;
t181 = t140 * t134;
t125 = atan2(t181, t135);
t122 = cos(t125);
t121 = sin(t125);
t165 = t121 * t181;
t107 = t122 * t135 + t165;
t104 = 0.1e1 / t107;
t105 = 0.1e1 / t107 ^ 2;
t123 = 0.1e1 / t126;
t200 = t123 - 0.1e1;
t138 = t142 ^ 2;
t186 = t130 * t138;
t101 = t105 * t186 + 0.1e1;
t175 = qJD(3) * t135;
t174 = qJD(3) * t140;
t162 = t132 * t174;
t163 = t134 * t176;
t98 = ((t135 * t174 + t163) * t131 + t130 * t162) * t123;
t157 = -t98 + t174;
t158 = t140 * t98 - qJD(3);
t188 = t122 * t134;
t93 = t158 * t188 + (t157 * t135 + t163) * t121;
t197 = t104 * t105 * t93;
t199 = 0.1e1 / t101 ^ 2 * (-t186 * t197 + (t134 * t138 * t175 - t153) * t105);
t99 = 0.1e1 / t101;
t198 = t105 * t99;
t191 = t114 * t139;
t189 = t121 * t135;
t183 = t134 * t142;
t172 = -0.2e1 * t199;
t171 = -0.2e1 * t197;
t170 = 0.2e1 * t202;
t169 = t104 * t199;
t168 = t99 * t175;
t166 = t103 * t116 * t119;
t164 = t123 * t130 * t131;
t159 = t131 * t203;
t154 = t140 * t164;
t152 = t160 * t142;
t151 = t141 * t190 + t191;
t118 = -t135 * t179 + t182;
t117 = t135 * t180 + t178;
t110 = 0.1e1 / t112;
t109 = t123 * t201;
t97 = (-t200 * t134 * t121 + t122 * t154) * t142;
t95 = t140 * t189 - t188 + (t122 * t181 - t189) * t109;
t94 = t201 * t203 + (qJD(1) * t152 + 0.2e1 * t140 * t150) * t123;
t1 = [t159 * t183 + (qJD(3) * t152 - t177 * t184) * t123, 0, t94, 0, 0, 0; (t104 * t168 + (-0.2e1 * t169 + (-qJD(1) * t97 - t93) * t198) * t134) * t140 + ((t97 * t168 + (t97 * t172 + ((0.2e1 * t134 * t194 - t98 * t154 - t200 * t175) * t121 + (t159 * t185 + t134 * t98 + (t129 * t162 + (-t98 + 0.2e1 * t174) * t134) * t123) * t122) * t99 * t142) * t134) * t105 + (t97 * t171 + (t104 + ((-t137 + t138) * t122 * t164 + t200 * t165) * t105) * qJD(1)) * t134 * t99) * t142, 0 (t104 * t99 * t177 + (0.2e1 * t169 + (qJD(3) * t95 + t93) * t198) * t142) * t135 + (((qJD(3) * t104 + t95 * t171) * t142 + (-t95 * t177 + ((t109 * t176 + t140 * t94) * t122 + ((-t109 * t140 + 0.1e1) * t98 + (t109 - t140) * qJD(3)) * t121) * t183) * t105) * t99 + (t95 * t172 + ((-t94 + t176) * t121 + (t157 * t109 + t158) * t122) * t99 * t135) * t105 * t142) * t134, 0, 0, 0; (-t114 * t117 + t118 * t190) * t170 + (-0.2e1 * t118 * t166 + t156 * t114 * t179 + (-t134 * t174 + t155 * t142) * t191 + (-t118 * t102 + t117 * t103 + (t155 * t178 - (qJD(3) * t134 * t141 + t156 * t139) * t140) * t119) * t115) * t110, 0, -0.2e1 * t151 * t183 * t202 + (t151 * t135 * t173 + (-t151 * t177 + ((qJD(4) * t114 + 0.2e1 * t166) * t141 + (t102 * t141 + (-qJD(4) * t119 + t103) * t139) * t115) * t142) * t134) * t110 (t114 * t120 + t193) * t170 + (-0.2e1 * t167 + (-t115 * t120 + t114 - 0.2e1 * t192) * t103) * t110, 0, 0;];
JaD_rot  = t1;
