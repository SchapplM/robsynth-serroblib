% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP7_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:19
% EndTime: 2019-02-26 21:11:20
% DurationCPUTime: 0.79s
% Computational Cost: add. (2887->96), mult. (2734->202), div. (498->12), fcn. (3199->9), ass. (0->96)
t161 = pkin(10) + qJ(3);
t157 = sin(t161);
t153 = t157 ^ 2;
t158 = cos(t161);
t155 = 0.1e1 / t158 ^ 2;
t211 = t153 * t155;
t166 = sin(qJ(1));
t229 = 0.2e1 * t166;
t228 = t157 * t211;
t163 = t166 ^ 2;
t149 = t163 * t211 + 0.1e1;
t147 = 0.1e1 / t149;
t154 = 0.1e1 / t158;
t167 = cos(qJ(1));
t201 = qJD(1) * t167;
t189 = t157 * t201;
t199 = qJD(3) * t166;
t120 = (-(-t158 * t199 - t189) * t154 + t199 * t211) * t147;
t227 = t120 - t199;
t165 = qJ(4) + qJ(5);
t159 = sin(t165);
t204 = t166 * t159;
t160 = cos(t165);
t206 = t160 * t167;
t142 = t158 * t206 + t204;
t205 = t166 * t157;
t145 = atan2(-t205, -t158);
t144 = cos(t145);
t143 = sin(t145);
t192 = t143 * t205;
t129 = -t144 * t158 - t192;
t126 = 0.1e1 / t129;
t136 = 0.1e1 / t142;
t127 = 0.1e1 / t129 ^ 2;
t137 = 0.1e1 / t142 ^ 2;
t226 = -0.2e1 * t157;
t225 = t147 - 0.1e1;
t213 = t144 * t157;
t115 = (-t120 * t166 + qJD(3)) * t213 + (t227 * t158 - t189) * t143;
t224 = t115 * t126 * t127;
t162 = qJD(4) + qJD(5);
t177 = t158 * t204 + t206;
t198 = qJD(3) * t167;
t188 = t157 * t198;
t121 = t177 * qJD(1) - t142 * t162 + t159 * t188;
t203 = t166 * t160;
t207 = t159 * t167;
t141 = t158 * t207 - t203;
t135 = t141 ^ 2;
t134 = t135 * t137 + 0.1e1;
t216 = t137 * t141;
t182 = -qJD(1) * t158 + t162;
t183 = t158 * t162 - qJD(1);
t122 = -t183 * t207 + (t182 * t166 - t188) * t160;
t221 = t122 * t136 * t137;
t223 = (-t121 * t216 - t135 * t221) / t134 ^ 2;
t222 = t120 * t157;
t220 = t127 * t157;
t219 = t127 * t167;
t209 = t154 * t157;
t176 = qJD(3) * (t154 * t228 + t209);
t180 = t153 * t166 * t201;
t218 = (t155 * t180 + t163 * t176) / t149 ^ 2;
t217 = t136 * t159;
t215 = t141 * t160;
t214 = t143 * t166;
t212 = t153 * t154;
t164 = t167 ^ 2;
t210 = t153 * t164;
t208 = t157 * t167;
t202 = qJD(1) * t166;
t200 = qJD(3) * t158;
t125 = t127 * t210 + 0.1e1;
t197 = 0.2e1 * (-t210 * t224 + (t157 * t164 * t200 - t180) * t127) / t125 ^ 2;
t196 = 0.2e1 * t224;
t195 = -0.2e1 * t223;
t194 = t141 * t221;
t193 = t127 * t208;
t191 = t147 * t212;
t187 = 0.1e1 + t211;
t186 = t157 * t197;
t185 = t218 * t226;
t184 = t218 * t229;
t181 = t166 * t191;
t179 = t187 * t167;
t178 = t137 * t215 - t217;
t175 = t157 * t199 + t182 * t167;
t140 = -t158 * t203 + t207;
t132 = 0.1e1 / t134;
t131 = t187 * t166 * t147;
t123 = 0.1e1 / t125;
t119 = (t225 * t157 * t143 - t144 * t181) * t167;
t118 = -t158 * t214 + t213 + (t143 * t158 - t144 * t205) * t131;
t116 = -t187 * t184 + (qJD(1) * t179 + t176 * t229) * t147;
t113 = t195 + 0.2e1 * (-t121 * t132 * t137 + (-t132 * t221 - t137 * t223) * t141) * t141;
t1 = [t154 * t167 * t185 + (qJD(3) * t179 - t202 * t209) * t147, 0, t116, 0, 0, 0; (t126 * t186 + (-t126 * t200 + (qJD(1) * t119 + t115) * t220) * t123) * t166 + (t127 * t186 * t119 + (-((t120 * t181 + t225 * t200 + t185) * t143 + (t184 * t212 - t222 + (t222 + (t226 - t228) * t199) * t147) * t144) * t193 + (-t127 * t200 + t157 * t196) * t119 + (-t126 + ((-t163 + t164) * t144 * t191 + t225 * t192) * t127) * t157 * qJD(1)) * t123) * t167, 0 (t118 * t220 - t126 * t158) * t167 * t197 + ((-t126 * t202 + (-qJD(3) * t118 - t115) * t219) * t158 + (-t126 * t198 - (-t116 * t144 * t166 - t227 * t143 + (-qJD(3) * t143 + t120 * t214 - t144 * t201) * t131) * t193 + (t127 * t202 + t167 * t196) * t118 - ((t116 - t201) * t143 + ((-t131 * t166 + 0.1e1) * qJD(3) + (t131 - t166) * t120) * t144) * t158 * t219) * t157) * t123, 0, 0, 0; 0.2e1 * (t136 * t177 + t140 * t216) * t223 + (0.2e1 * t140 * t194 - t183 * t136 * t203 + t175 * t217 + (-t183 * t141 * t204 + t140 * t121 + t122 * t177 - t175 * t215) * t137) * t132, 0, t178 * t195 * t208 + (t178 * t158 * t198 + (-t178 * t202 + ((-t136 * t162 - 0.2e1 * t194) * t160 + (-t121 * t160 + (-t141 * t162 + t122) * t159) * t137) * t167) * t157) * t132, t113, t113, 0;];
JaD_rot  = t1;
