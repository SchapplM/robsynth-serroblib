% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRR9
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRR9_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_jacobiaD_rot_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:19:14
% EndTime: 2019-02-26 21:19:15
% DurationCPUTime: 0.73s
% Computational Cost: add. (2159->94), mult. (2949->205), div. (516->12), fcn. (3430->9), ass. (0->97)
t157 = sin(qJ(3));
t151 = 0.1e1 / t157 ^ 2;
t159 = cos(qJ(3));
t155 = t159 ^ 2;
t204 = t151 * t155;
t160 = cos(qJ(1));
t222 = 0.2e1 * t160;
t221 = t159 * t204;
t197 = t160 * t159;
t141 = atan2(-t197, t157);
t139 = sin(t141);
t140 = cos(t141);
t128 = -t139 * t197 + t140 * t157;
t125 = 0.1e1 / t128;
t149 = qJ(4) + qJ(5) + qJ(6);
t146 = sin(t149);
t147 = cos(t149);
t158 = sin(qJ(1));
t200 = t158 * t147;
t136 = t146 * t160 + t157 * t200;
t132 = 0.1e1 / t136;
t150 = 0.1e1 / t157;
t126 = 0.1e1 / t128 ^ 2;
t133 = 0.1e1 / t136 ^ 2;
t156 = t160 ^ 2;
t144 = t156 * t204 + 0.1e1;
t142 = 0.1e1 / t144;
t220 = t142 - 0.1e1;
t153 = t158 ^ 2;
t203 = t153 * t155;
t124 = t126 * t203 + 0.1e1;
t194 = qJD(1) * t160;
t175 = t155 * t158 * t194;
t192 = qJD(3) * t159;
t195 = qJD(1) * t159;
t184 = t158 * t195;
t191 = qJD(3) * t160;
t118 = ((t157 * t191 + t184) * t150 + t191 * t204) * t142;
t206 = t140 * t159;
t112 = (-t118 * t160 + qJD(3)) * t206 + (t184 + (-t118 + t191) * t157) * t139;
t217 = t112 * t125 * t126;
t219 = (-t203 * t217 + (-t153 * t157 * t192 + t175) * t126) / t124 ^ 2;
t148 = qJD(4) + qJD(5) + qJD(6);
t178 = qJD(1) * t157 + t148;
t168 = t158 * t192 + t178 * t160;
t179 = t148 * t157 + qJD(1);
t172 = t147 * t179;
t116 = t168 * t146 + t158 * t172;
t198 = t160 * t147;
t201 = t158 * t146;
t135 = t157 * t201 - t198;
t131 = t135 ^ 2;
t121 = t131 * t133 + 0.1e1;
t208 = t133 * t135;
t173 = t146 * t179;
t117 = t168 * t147 - t158 * t173;
t216 = t117 * t132 * t133;
t218 = (t116 * t208 - t131 * t216) / t121 ^ 2;
t215 = t118 * t139;
t214 = t118 * t159;
t213 = t126 * t158;
t212 = t126 * t159;
t170 = qJD(3) * (-t159 - t221) * t150;
t211 = (-t151 * t175 + t156 * t170) / t144 ^ 2;
t182 = 0.1e1 + t204;
t130 = t182 * t160 * t142;
t210 = t130 * t160;
t209 = t132 * t146;
t207 = t135 * t147;
t205 = t150 * t155;
t202 = t157 * t160;
t199 = t158 * t159;
t196 = qJD(1) * t158;
t193 = qJD(3) * t157;
t190 = 0.2e1 * t218;
t189 = -0.2e1 * t217;
t188 = t159 * t219;
t187 = t126 * t199;
t186 = t159 * t211;
t185 = t142 * t205;
t183 = t159 * t194;
t181 = 0.2e1 * t135 * t216;
t180 = t211 * t222;
t177 = t160 * t185;
t176 = t220 * t159 * t139;
t174 = t182 * t158;
t171 = t133 * t207 - t209;
t169 = -t178 * t158 + t159 * t191;
t138 = t157 * t198 - t201;
t137 = t146 * t202 + t200;
t122 = 0.1e1 / t124;
t119 = 0.1e1 / t121;
t115 = (-t140 * t177 - t176) * t158;
t114 = t139 * t202 + t206 + (-t139 * t157 - t140 * t197) * t130;
t113 = -t182 * t180 + (-qJD(1) * t174 + t170 * t222) * t142;
t109 = -0.2e1 * t218 + 0.2e1 * (t116 * t119 * t133 + (-t119 * t216 - t133 * t218) * t135) * t135;
t1 = [-0.2e1 * t158 * t150 * t186 + (-qJD(3) * t174 + t150 * t183) * t142, 0, t113, 0, 0, 0; (0.2e1 * t125 * t188 + (t125 * t193 + (qJD(1) * t115 + t112) * t212) * t122) * t160 + (-0.2e1 * t126 * t188 * t115 + (((t118 * t177 + t220 * t193 + 0.2e1 * t186) * t139 + (t180 * t205 + t214 + (-t214 + (0.2e1 * t159 + t221) * t191) * t142) * t140) * t187 + (-t126 * t193 + t159 * t189) * t115 + (t125 + ((t153 - t156) * t140 * t185 - t160 * t176) * t126) * t195) * t122) * t158, 0, 0.2e1 * (-t114 * t212 - t125 * t157) * t158 * t219 + ((t125 * t194 + (-qJD(3) * t114 - t112) * t213) * t157 + (t158 * qJD(3) * t125 + (-t113 * t140 * t160 + t139 * t191 + t210 * t215 - t215 + (-qJD(3) * t139 + t140 * t196) * t130) * t187 + (t126 * t194 + t158 * t189) * t114 + ((-t113 - t196) * t139 + ((-0.1e1 + t210) * qJD(3) + (-t130 + t160) * t118) * t140) * t157 * t213) * t159) * t122, 0, 0, 0; (-t132 * t137 + t138 * t208) * t190 + (t138 * t181 + t160 * t132 * t172 + t169 * t209 + (t160 * t135 * t173 - t138 * t116 - t137 * t117 - t169 * t207) * t133) * t119, 0, t171 * t190 * t199 + (-t171 * t183 + (t171 * t193 + ((t132 * t148 + t181) * t147 + (-t116 * t147 + (t135 * t148 - t117) * t146) * t133) * t159) * t158) * t119, t109, t109, t109;];
JaD_rot  = t1;
