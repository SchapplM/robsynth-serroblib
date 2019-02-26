% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRR1_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_jacobiaD_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:43
% EndTime: 2019-02-26 20:34:44
% DurationCPUTime: 0.75s
% Computational Cost: add. (7044->97), mult. (3810->206), div. (753->12), fcn. (4455->9), ass. (0->96)
t180 = qJ(1) + pkin(10);
t176 = sin(t180);
t243 = 0.2e1 * t176;
t174 = t176 ^ 2;
t178 = pkin(11) + qJ(4) + qJ(5);
t172 = sin(t178);
t168 = t172 ^ 2;
t173 = cos(t178);
t170 = 0.1e1 / t173 ^ 2;
t229 = t168 * t170;
t163 = t174 * t229 + 0.1e1;
t161 = 0.1e1 / t163;
t169 = 0.1e1 / t173;
t177 = cos(t180);
t214 = qJD(1) * t177;
t203 = t172 * t214;
t179 = qJD(4) + qJD(5);
t221 = t176 * t179;
t207 = t170 * t221;
t135 = (-(-t173 * t221 - t203) * t169 + t168 * t207) * t161;
t242 = t135 - t221;
t182 = cos(qJ(6));
t217 = t177 * t182;
t181 = sin(qJ(6));
t220 = t176 * t181;
t157 = t173 * t217 + t220;
t198 = qJD(6) * t173 - qJD(1);
t224 = t172 * t179;
t241 = t181 * t198 + t182 * t224;
t222 = t176 * t172;
t160 = atan2(-t222, -t173);
t159 = cos(t160);
t158 = sin(t160);
t208 = t158 * t222;
t145 = -t159 * t173 - t208;
t142 = 0.1e1 / t145;
t151 = 0.1e1 / t157;
t143 = 0.1e1 / t145 ^ 2;
t152 = 0.1e1 / t157 ^ 2;
t240 = t161 - 0.1e1;
t231 = t159 * t172;
t130 = (-t135 * t176 + t179) * t231 + (t242 * t173 - t203) * t158;
t239 = t130 * t142 * t143;
t191 = t173 * t220 + t217;
t206 = t181 * t224;
t139 = t191 * qJD(1) - t157 * qJD(6) + t177 * t206;
t218 = t177 * t181;
t219 = t176 * t182;
t156 = t173 * t218 - t219;
t150 = t156 ^ 2;
t149 = t150 * t152 + 0.1e1;
t233 = t152 * t156;
t197 = -qJD(1) * t173 + qJD(6);
t193 = t197 * t182;
t140 = t176 * t193 - t241 * t177;
t237 = t140 * t151 * t152;
t238 = (-t139 * t233 - t150 * t237) / t149 ^ 2;
t167 = t172 * t168;
t226 = t169 * t172;
t190 = (t167 * t169 * t170 + t226) * t179;
t227 = t168 * t176;
t195 = t214 * t227;
t236 = (t170 * t195 + t174 * t190) / t163 ^ 2;
t235 = t143 * t172;
t234 = t143 * t177;
t232 = t158 * t176;
t230 = t168 * t169;
t175 = t177 ^ 2;
t228 = t168 * t175;
t225 = t172 * t177;
t223 = t173 * t179;
t216 = t179 * t142;
t215 = qJD(1) * t176;
t138 = t143 * t228 + 0.1e1;
t213 = 0.2e1 * (-t228 * t239 + (t172 * t175 * t223 - t195) * t143) / t138 ^ 2;
t212 = 0.2e1 * t239;
t211 = -0.2e1 * t238;
t210 = t143 * t225;
t209 = t156 * t237;
t202 = 0.1e1 + t229;
t201 = t172 * t213;
t200 = -0.2e1 * t172 * t236;
t199 = t236 * t243;
t196 = t159 * t161 * t230;
t194 = t202 * t177;
t192 = -t151 * t181 + t182 * t233;
t155 = -t173 * t219 + t218;
t147 = 0.1e1 / t149;
t146 = t202 * t176 * t161;
t136 = 0.1e1 / t138;
t134 = (t158 * t172 * t240 - t176 * t196) * t177;
t132 = -t173 * t232 + t231 + (t158 * t173 - t159 * t222) * t146;
t131 = -t202 * t199 + (qJD(1) * t194 + t190 * t243) * t161;
t128 = t192 * t211 * t225 + (t192 * t177 * t223 + (-t192 * t215 + ((-qJD(6) * t151 - 0.2e1 * t209) * t182 + (-t139 * t182 + (-qJD(6) * t156 + t140) * t181) * t152) * t177) * t172) * t147;
t127 = (t132 * t235 - t142 * t173) * t177 * t213 + ((-t142 * t215 + (-t132 * t179 - t130) * t234) * t173 + (-t177 * t216 - (-t131 * t159 * t176 - t242 * t158 + (t135 * t232 - t158 * t179 - t159 * t214) * t146) * t210 + (t143 * t215 + t177 * t212) * t132 - ((t131 - t214) * t158 + ((-t146 * t176 + 0.1e1) * t179 + (t146 - t176) * t135) * t159) * t173 * t234) * t172) * t136;
t1 = [t177 * t169 * t200 + (t179 * t194 - t215 * t226) * t161, 0, 0, t131, t131, 0; (t142 * t201 + (-t173 * t216 + (qJD(1) * t134 + t130) * t235) * t136) * t176 + (t143 * t201 * t134 + (-((t200 - t223 + (t135 * t169 * t227 + t223) * t161) * t158 + (t199 * t230 - t135 * t172 + (-t167 * t207 + (t135 - 0.2e1 * t221) * t172) * t161) * t159) * t210 + (-t143 * t223 + t172 * t212) * t134 + (-t142 + ((-t174 + t175) * t196 + t240 * t208) * t143) * t172 * qJD(1)) * t136) * t177, 0, 0, t127, t127, 0; 0.2e1 * (t151 * t191 + t155 * t233) * t238 + (0.2e1 * t155 * t209 + (t155 * t139 + t191 * t140 + (-t241 * t176 - t177 * t193) * t156) * t152 + (t197 * t218 + (-t182 * t198 + t206) * t176) * t151) * t147, 0, 0, t128, t128, t211 + 0.2e1 * (-t139 * t152 * t147 + (-t147 * t237 - t152 * t238) * t156) * t156;];
JaD_rot  = t1;
