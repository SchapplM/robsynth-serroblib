% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPR6_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR6_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR6_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:47
% EndTime: 2019-02-26 21:03:48
% DurationCPUTime: 0.78s
% Computational Cost: add. (3326->96), mult. (2734->203), div. (498->12), fcn. (3199->9), ass. (0->95)
t162 = pkin(10) + qJ(3);
t159 = sin(t162);
t155 = t159 ^ 2;
t160 = cos(t162);
t157 = 0.1e1 / t160 ^ 2;
t210 = t155 * t157;
t166 = sin(qJ(1));
t228 = 0.2e1 * t166;
t227 = t159 * t210;
t164 = t166 ^ 2;
t149 = t164 * t210 + 0.1e1;
t147 = 0.1e1 / t149;
t156 = 0.1e1 / t160;
t167 = cos(qJ(1));
t201 = qJD(1) * t167;
t189 = t159 * t201;
t199 = qJD(3) * t166;
t122 = (-(-t160 * t199 - t189) * t156 + t199 * t210) * t147;
t226 = t122 - t199;
t161 = qJ(4) + pkin(11) + qJ(6);
t153 = cos(t161);
t152 = sin(t161);
t205 = t166 * t152;
t206 = t160 * t167;
t142 = t153 * t206 + t205;
t203 = t166 * t159;
t146 = atan2(-t203, -t160);
t145 = cos(t146);
t144 = sin(t146);
t192 = t144 * t203;
t132 = -t145 * t160 - t192;
t129 = 0.1e1 / t132;
t136 = 0.1e1 / t142;
t130 = 0.1e1 / t132 ^ 2;
t137 = 0.1e1 / t142 ^ 2;
t225 = -0.2e1 * t159;
t224 = t147 - 0.1e1;
t213 = t145 * t159;
t115 = (-t122 * t166 + qJD(3)) * t213 + (t226 * t160 - t189) * t144;
t223 = t115 * t129 * t130;
t163 = qJD(4) + qJD(6);
t177 = t153 * t167 + t160 * t205;
t198 = qJD(3) * t167;
t188 = t159 * t198;
t120 = t177 * qJD(1) - t142 * t163 + t152 * t188;
t204 = t166 * t153;
t141 = t152 * t206 - t204;
t135 = t141 ^ 2;
t128 = t135 * t137 + 0.1e1;
t216 = t137 * t141;
t182 = -qJD(1) * t160 + t163;
t183 = t160 * t163 - qJD(1);
t212 = t152 * t167;
t121 = -t183 * t212 + (t182 * t166 - t188) * t153;
t221 = t121 * t136 * t137;
t222 = (-t120 * t216 - t135 * t221) / t128 ^ 2;
t220 = t122 * t159;
t219 = t130 * t159;
t208 = t156 * t159;
t176 = qJD(3) * (t156 * t227 + t208);
t180 = t155 * t166 * t201;
t218 = (t157 * t180 + t164 * t176) / t149 ^ 2;
t217 = t136 * t152;
t215 = t141 * t153;
t214 = t144 * t166;
t211 = t155 * t156;
t165 = t167 ^ 2;
t209 = t155 * t165;
t207 = t159 * t167;
t202 = qJD(1) * t166;
t200 = qJD(3) * t160;
t125 = t130 * t209 + 0.1e1;
t197 = 0.2e1 * (-t209 * t223 + (t159 * t165 * t200 - t180) * t130) / t125 ^ 2;
t196 = 0.2e1 * t223;
t195 = -0.2e1 * t222;
t194 = t141 * t221;
t193 = t130 * t207;
t191 = t147 * t211;
t187 = 0.1e1 + t210;
t186 = t159 * t197;
t185 = t218 * t225;
t184 = t218 * t228;
t181 = t166 * t191;
t179 = t187 * t167;
t178 = t137 * t215 - t217;
t175 = t159 * t199 + t182 * t167;
t140 = -t160 * t204 + t212;
t134 = t187 * t166 * t147;
t126 = 0.1e1 / t128;
t123 = 0.1e1 / t125;
t119 = (t224 * t159 * t144 - t145 * t181) * t167;
t118 = -t160 * t214 + t213 + (t144 * t160 - t145 * t203) * t134;
t116 = -t187 * t184 + (qJD(1) * t179 + t176 * t228) * t147;
t113 = t195 + 0.2e1 * (-t120 * t126 * t137 + (-t126 * t221 - t137 * t222) * t141) * t141;
t1 = [t156 * t167 * t185 + (qJD(3) * t179 - t202 * t208) * t147, 0, t116, 0, 0, 0; (t129 * t186 + (-t129 * t200 + (qJD(1) * t119 + t115) * t219) * t123) * t166 + (t130 * t186 * t119 + (-((t122 * t181 + t224 * t200 + t185) * t144 + (t184 * t211 - t220 + (t220 + (t225 - t227) * t199) * t147) * t145) * t193 + (-t130 * t200 + t159 * t196) * t119 + (-t129 + ((-t164 + t165) * t145 * t191 + t224 * t192) * t130) * t159 * qJD(1)) * t123) * t167, 0 (t118 * t219 - t129 * t160) * t167 * t197 + ((-t129 * t202 + (-qJD(3) * t118 - t115) * t167 * t130) * t160 + (-t129 * t198 - (-t116 * t145 * t166 - t226 * t144 + (-qJD(3) * t144 + t122 * t214 - t145 * t201) * t134) * t193 + (t130 * t202 + t167 * t196) * t118 - ((t116 - t201) * t144 + ((-t134 * t166 + 0.1e1) * qJD(3) + (t134 - t166) * t122) * t145) * t130 * t206) * t159) * t123, 0, 0, 0; 0.2e1 * (t136 * t177 + t140 * t216) * t222 + (0.2e1 * t140 * t194 - t183 * t136 * t204 + t175 * t217 + (-t183 * t141 * t205 + t140 * t120 + t121 * t177 - t175 * t215) * t137) * t126, 0, t178 * t195 * t207 + (t178 * t160 * t198 + (-t178 * t202 + ((-t136 * t163 - 0.2e1 * t194) * t153 + (-t120 * t153 + (-t141 * t163 + t121) * t152) * t137) * t167) * t159) * t126, t113, 0, t113;];
JaD_rot  = t1;
