% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRRRR4_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:47
% EndTime: 2019-02-26 22:48:48
% DurationCPUTime: 0.70s
% Computational Cost: add. (3336->96), mult. (3164->203), div. (534->12), fcn. (3661->9), ass. (0->96)
t164 = sin(qJ(2));
t158 = t164 ^ 2;
t166 = cos(qJ(2));
t161 = 0.1e1 / t166 ^ 2;
t212 = t158 * t161;
t165 = sin(qJ(1));
t230 = 0.2e1 * t165;
t229 = t164 * t212;
t156 = qJ(3) + qJ(4) + qJ(5) + qJ(6);
t154 = cos(t156);
t167 = cos(qJ(1));
t203 = t167 * t154;
t153 = sin(t156);
t208 = t165 * t153;
t143 = t166 * t203 + t208;
t206 = t165 * t164;
t148 = atan2(-t206, -t166);
t147 = cos(t148);
t146 = sin(t148);
t192 = t146 * t206;
t133 = -t147 * t166 - t192;
t130 = 0.1e1 / t133;
t137 = 0.1e1 / t143;
t160 = 0.1e1 / t166;
t131 = 0.1e1 / t133 ^ 2;
t138 = 0.1e1 / t143 ^ 2;
t228 = -0.2e1 * t164;
t159 = t165 ^ 2;
t152 = t159 * t212 + 0.1e1;
t149 = 0.1e1 / t152;
t227 = t149 - 0.1e1;
t155 = qJD(3) + qJD(4) + qJD(5) + qJD(6);
t205 = t165 * t166;
t177 = t153 * t205 + t203;
t198 = qJD(2) * t167;
t188 = t164 * t198;
t121 = t177 * qJD(1) - t143 * t155 + t153 * t188;
t204 = t167 * t153;
t207 = t165 * t154;
t142 = t166 * t204 - t207;
t136 = t142 ^ 2;
t125 = t136 * t138 + 0.1e1;
t216 = t138 * t142;
t182 = -qJD(1) * t166 + t155;
t183 = t155 * t166 - qJD(1);
t122 = -t183 * t204 + (t182 * t165 - t188) * t154;
t224 = t122 * t137 * t138;
t226 = (-t121 * t216 - t136 * t224) / t125 ^ 2;
t201 = qJD(1) * t167;
t189 = t164 * t201;
t199 = qJD(2) * t166;
t200 = qJD(2) * t165;
t126 = (-(-t165 * t199 - t189) * t160 + t200 * t212) * t149;
t214 = t147 * t164;
t117 = (-t126 * t165 + qJD(2)) * t214 + (-t189 + (t126 - t200) * t166) * t146;
t225 = t117 * t130 * t131;
t223 = t126 * t146;
t222 = t126 * t164;
t221 = t131 * t164;
t220 = t131 * t167;
t210 = t160 * t164;
t176 = qJD(2) * (t160 * t229 + t210);
t180 = t158 * t165 * t201;
t219 = (t159 * t176 + t161 * t180) / t152 ^ 2;
t187 = 0.1e1 + t212;
t135 = t187 * t165 * t149;
t218 = t135 * t165;
t217 = t137 * t153;
t215 = t142 * t154;
t213 = t158 * t160;
t163 = t167 ^ 2;
t211 = t158 * t163;
t209 = t164 * t167;
t202 = qJD(1) * t165;
t129 = t131 * t211 + 0.1e1;
t197 = 0.2e1 * (-t211 * t225 + (t163 * t164 * t199 - t180) * t131) / t129 ^ 2;
t196 = -0.2e1 * t226;
t195 = 0.2e1 * t225;
t194 = t131 * t209;
t193 = t142 * t224;
t191 = t149 * t213;
t186 = t164 * t197;
t185 = t219 * t228;
t184 = t219 * t230;
t181 = t165 * t191;
t179 = t187 * t167;
t178 = t138 * t215 - t217;
t175 = t164 * t200 + t182 * t167;
t141 = -t154 * t205 + t204;
t127 = 0.1e1 / t129;
t123 = 0.1e1 / t125;
t120 = (t227 * t164 * t146 - t147 * t181) * t167;
t119 = -t146 * t205 + t214 + (t146 * t166 - t147 * t206) * t135;
t118 = -t187 * t184 + (qJD(1) * t179 + t176 * t230) * t149;
t114 = t196 + 0.2e1 * (-t121 * t138 * t123 + (-t123 * t224 - t138 * t226) * t142) * t142;
t1 = [t167 * t160 * t185 + (qJD(2) * t179 - t202 * t210) * t149, t118, 0, 0, 0, 0; (t130 * t186 + (-t130 * t199 + (qJD(1) * t120 + t117) * t221) * t127) * t165 + (t131 * t186 * t120 + (-((t126 * t181 + t227 * t199 + t185) * t146 + (t184 * t213 - t222 + (t222 + (t228 - t229) * t200) * t149) * t147) * t194 + (-t131 * t199 + t164 * t195) * t120 + (-t130 + ((-t159 + t163) * t147 * t191 + t227 * t192) * t131) * t164 * qJD(1)) * t127) * t167 (t119 * t221 - t130 * t166) * t167 * t197 + ((-t130 * t202 + (-qJD(2) * t119 - t117) * t220) * t166 + (-t130 * t198 - (-t118 * t147 * t165 + t146 * t200 + t218 * t223 - t223 + (-qJD(2) * t146 - t147 * t201) * t135) * t194 + (t131 * t202 + t167 * t195) * t119 - ((t118 - t201) * t146 + ((0.1e1 - t218) * qJD(2) + (t135 - t165) * t126) * t147) * t166 * t220) * t164) * t127, 0, 0, 0, 0; 0.2e1 * (t137 * t177 + t141 * t216) * t226 + (0.2e1 * t141 * t193 - t183 * t137 * t207 + t175 * t217 + (-t183 * t142 * t208 + t141 * t121 + t122 * t177 - t175 * t215) * t138) * t123, t178 * t196 * t209 + (t178 * t166 * t198 + (-t178 * t202 + ((-t137 * t155 - 0.2e1 * t193) * t154 + (-t121 * t154 + (-t142 * t155 + t122) * t153) * t138) * t167) * t164) * t123, t114, t114, t114, t114;];
JaD_rot  = t1;
