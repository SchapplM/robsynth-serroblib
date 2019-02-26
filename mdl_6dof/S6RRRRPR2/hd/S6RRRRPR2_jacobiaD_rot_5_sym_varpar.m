% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:16
% EndTime: 2019-02-26 22:31:17
% DurationCPUTime: 0.73s
% Computational Cost: add. (7384->87), mult. (4637->191), div. (990->12), fcn. (5438->9), ass. (0->91)
t157 = sin(qJ(1));
t212 = 0.2e1 * t157;
t152 = qJ(2) + qJ(3) + qJ(4);
t150 = cos(t152);
t155 = sin(pkin(11));
t158 = cos(qJ(1));
t187 = t158 * t155;
t156 = cos(pkin(11));
t188 = t157 * t156;
t138 = t150 * t187 - t188;
t130 = t138 ^ 2;
t186 = t158 * t156;
t189 = t157 * t155;
t139 = t150 * t186 + t189;
t132 = 0.1e1 / t139 ^ 2;
t127 = t130 * t132 + 0.1e1;
t136 = -t150 * t189 - t186;
t149 = sin(t152);
t151 = qJD(2) + qJD(3) + qJD(4);
t191 = t151 * t158;
t176 = t149 * t191;
t128 = t136 * qJD(1) - t155 * t176;
t202 = t132 * t138;
t137 = -t150 * t188 + t187;
t129 = t137 * qJD(1) - t156 * t176;
t131 = 0.1e1 / t139;
t204 = t129 * t131 * t132;
t211 = (t128 * t202 - t130 * t204) / t127 ^ 2;
t153 = t157 ^ 2;
t145 = t149 ^ 2;
t147 = 0.1e1 / t150 ^ 2;
t198 = t145 * t147;
t143 = t153 * t198 + 0.1e1;
t141 = 0.1e1 / t143;
t146 = 0.1e1 / t150;
t184 = qJD(1) * t158;
t175 = t149 * t184;
t192 = t151 * t157;
t178 = t147 * t192;
t115 = (-(-t150 * t192 - t175) * t146 + t145 * t178) * t141;
t210 = t115 - t192;
t190 = t157 * t149;
t140 = atan2(-t190, -t150);
t135 = cos(t140);
t134 = sin(t140);
t180 = t134 * t190;
t123 = -t135 * t150 - t180;
t120 = 0.1e1 / t123;
t121 = 0.1e1 / t123 ^ 2;
t209 = t141 - 0.1e1;
t200 = t135 * t149;
t110 = (-t115 * t157 + t151) * t200 + (t210 * t150 - t175) * t134;
t208 = t110 * t120 * t121;
t144 = t149 * t145;
t195 = t146 * t149;
t166 = t151 * (t144 * t146 * t147 + t195);
t196 = t145 * t157;
t169 = t184 * t196;
t207 = (t147 * t169 + t153 * t166) / t143 ^ 2;
t206 = t121 * t149;
t205 = t121 * t158;
t203 = t131 * t155;
t201 = t134 * t157;
t199 = t145 * t146;
t154 = t158 ^ 2;
t197 = t145 * t154;
t194 = t149 * t158;
t193 = t150 * t151;
t185 = qJD(1) * t157;
t118 = t121 * t197 + 0.1e1;
t183 = 0.2e1 * (-t197 * t208 + (t149 * t154 * t193 - t169) * t121) / t118 ^ 2;
t182 = 0.2e1 * t208;
t181 = t121 * t194;
t179 = t138 * t204;
t177 = t151 * t190;
t174 = 0.1e1 + t198;
t173 = t149 * t183;
t172 = -0.2e1 * t149 * t207;
t171 = t207 * t212;
t170 = t135 * t141 * t199;
t168 = t174 * t158;
t167 = t156 * t202 - t203;
t125 = 0.1e1 / t127;
t124 = t174 * t157 * t141;
t116 = 0.1e1 / t118;
t114 = (t209 * t149 * t134 - t157 * t170) * t158;
t112 = -t150 * t201 + t200 + (t134 * t150 - t135 * t190) * t124;
t111 = -t174 * t171 + (qJD(1) * t168 + t166 * t212) * t141;
t108 = -0.2e1 * t167 * t194 * t211 + (t167 * t150 * t191 + (-0.2e1 * t179 * t186 + t185 * t203 + (t129 * t187 + (t128 * t158 - t138 * t185) * t156) * t132) * t149) * t125;
t107 = (t112 * t206 - t120 * t150) * t158 * t183 + ((-t120 * t185 + (-t112 * t151 - t110) * t205) * t150 + (-t120 * t191 - (-t111 * t135 * t157 - t210 * t134 + (t115 * t201 - t134 * t151 - t135 * t184) * t124) * t181 + (t121 * t185 + t158 * t182) * t112 - ((t111 - t184) * t134 + ((-t124 * t157 + 0.1e1) * t151 + (t124 - t157) * t115) * t135) * t150 * t205) * t149) * t116;
t1 = [t158 * t146 * t172 + (t151 * t168 - t185 * t195) * t141, t111, t111, t111, 0, 0; (t120 * t173 + (-t120 * t193 + (qJD(1) * t114 + t110) * t206) * t116) * t157 + (t121 * t173 * t114 + (-((t172 - t193 + (t115 * t146 * t196 + t193) * t141) * t134 + (t171 * t199 - t115 * t149 + (-t144 * t178 + (t115 - 0.2e1 * t192) * t149) * t141) * t135) * t181 + (-t121 * t193 + t149 * t182) * t114 + (-t120 + ((-t153 + t154) * t170 + t209 * t180) * t121) * t149 * qJD(1)) * t116) * t158, t107, t107, t107, 0, 0; 0.2e1 * (-t131 * t136 + t137 * t202) * t211 + ((-t138 * qJD(1) + t155 * t177) * t131 + 0.2e1 * t137 * t179 + (-t136 * t129 - (-t139 * qJD(1) + t156 * t177) * t138 - t137 * t128) * t132) * t125, t108, t108, t108, 0, 0;];
JaD_rot  = t1;
