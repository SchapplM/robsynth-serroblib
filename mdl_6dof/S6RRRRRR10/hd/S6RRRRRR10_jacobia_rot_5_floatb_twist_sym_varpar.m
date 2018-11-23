% Rotatorische Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR10
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
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d1,d2,d3,d4,d5,d6]';
%
% Output:
% Ja_rot [3x6]
%   Rotatorische Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 11:27
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Ja_rot = S6RRRRRR10_jacobia_rot_5_floatb_twist_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(14,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10_jacobia_rot_5_floatb_twist_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6RRRRRR10_jacobia_rot_5_floatb_twist_sym_varpar: pkin has to be [14x1] (double)');

%% Symbolic Calculation
% From jacobia_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 11:27:18
% EndTime: 2018-11-23 11:27:19
% DurationCPUTime: 1.07s
% Computational Cost: add. (13311->104), mult. (12892->191), div. (115->9), fcn. (12904->31), ass. (0->108)
t254 = pkin(6) + qJ(2);
t241 = cos(t254) / 0.2e1;
t255 = pkin(6) - qJ(2);
t247 = cos(t255);
t203 = t247 / 0.2e1 + t241;
t224 = sin(qJ(2));
t225 = sin(qJ(1));
t230 = cos(qJ(1));
t187 = -t230 * t203 + t225 * t224;
t238 = sin(t254) / 0.2e1;
t244 = sin(t255);
t197 = t238 - t244 / 0.2e1;
t229 = cos(qJ(2));
t188 = t230 * t197 + t225 * t229;
t252 = pkin(7) + qJ(3);
t237 = sin(t252) / 0.2e1;
t253 = pkin(7) - qJ(3);
t243 = sin(t253);
t195 = t237 - t243 / 0.2e1;
t240 = cos(t252) / 0.2e1;
t246 = cos(t253);
t200 = t240 - t246 / 0.2e1;
t228 = cos(qJ(3));
t218 = sin(pkin(6));
t256 = t218 * t230;
t173 = -t187 * t195 + t188 * t228 + t200 * t256;
t194 = t237 + t243 / 0.2e1;
t201 = t246 / 0.2e1 + t240;
t223 = sin(qJ(3));
t174 = t187 * t201 + t188 * t223 + t194 * t256;
t250 = pkin(8) + qJ(4);
t236 = sin(t250) / 0.2e1;
t251 = pkin(8) - qJ(4);
t242 = sin(t251);
t193 = t236 - t242 / 0.2e1;
t239 = cos(t250) / 0.2e1;
t245 = cos(t251);
t198 = t239 - t245 / 0.2e1;
t227 = cos(qJ(4));
t217 = sin(pkin(7));
t267 = cos(pkin(7));
t248 = t218 * t267;
t232 = t187 * t217 - t230 * t248;
t153 = t173 * t227 - t174 * t193 - t232 * t198;
t199 = t245 / 0.2e1 + t239;
t222 = sin(qJ(4));
t234 = t236 + t242 / 0.2e1;
t151 = t173 * t222 + t174 * t199 - t232 * t234;
t196 = t238 + t244 / 0.2e1;
t202 = t241 - t247 / 0.2e1;
t220 = cos(pkin(6));
t181 = t220 * t194 + t196 * t201 + t202 * t223;
t182 = t196 * t195 - t220 * t200 - t202 * t228;
t186 = -t196 * t217 + t220 * t267;
t163 = -t181 * t199 + t182 * t222 - t186 * t234;
t146 = atan2(-t151, t163);
t143 = sin(t146);
t144 = cos(t146);
t141 = -t143 * t151 + t144 * t163;
t140 = 0.1e1 / t141 ^ 2;
t190 = -t225 * t203 - t230 * t224;
t192 = t225 * t197 - t230 * t229;
t257 = t218 * t225;
t176 = t190 * t201 + t192 * t223 + t194 * t257;
t178 = -t190 * t195 + t192 * t228 + t200 * t257;
t231 = t190 * t217 - t225 * t248;
t155 = -t176 * t199 - t178 * t222 + t231 * t234;
t266 = t140 * t155;
t265 = t144 * t151;
t156 = t176 * t193 - t178 * t227 + t231 * t198;
t216 = sin(pkin(8));
t219 = cos(pkin(8));
t168 = -t176 * t216 - t231 * t219;
t221 = sin(qJ(5));
t226 = cos(qJ(5));
t150 = t156 * t226 + t168 * t221;
t148 = 0.1e1 / t150 ^ 2;
t149 = t156 * t221 - t168 * t226;
t264 = t148 * t149;
t162 = 0.1e1 / t163 ^ 2;
t263 = t151 * t162;
t262 = t155 ^ 2 * t140;
t259 = t178 * t216;
t258 = t192 * t217;
t249 = t149 ^ 2 * t148 + 0.1e1;
t235 = -t143 * t163 - t265;
t233 = t217 * t234;
t180 = t190 * t228 + t192 * t195;
t179 = -t190 * t223 + t192 * t201;
t169 = -t179 * t216 - t219 * t258;
t167 = -t174 * t216 - t219 * t232;
t166 = (t202 * t195 + t196 * t228) * t222 - (-t196 * t223 + t202 * t201) * t199 + t202 * t233;
t165 = t181 * t222 + t182 * t199;
t164 = t181 * t193 + t182 * t227 - t186 * t198;
t161 = 0.1e1 / t163;
t160 = t176 * t227 + t178 * t193;
t159 = t173 * t199 - t174 * t222;
t158 = t179 * t193 + t180 * t227 + t198 * t258;
t157 = (-t187 * t228 - t188 * t195) * t222 - (t187 * t223 - t188 * t201) * t199 - t188 * t233;
t147 = 0.1e1 / t150;
t145 = 0.1e1 / (t151 ^ 2 * t162 + 0.1e1);
t142 = 0.1e1 / t249;
t139 = 0.1e1 / t141;
t138 = 0.1e1 / (0.1e1 + t262);
t137 = (-t157 * t161 + t166 * t263) * t145;
t136 = (-t159 * t161 + t165 * t263) * t145;
t135 = (-t153 * t161 + t164 * t263) * t145;
t1 = [-t155 * t161 * t145, t137, t136, t135, 0, 0; (-t151 * t139 - (-t143 + (t161 * t265 + t143) * t145) * t262) * t138 ((-t179 * t199 + t180 * t222 + t192 * t233) * t139 - (t235 * t137 - t143 * t157 + t144 * t166) * t266) * t138 ((t176 * t222 - t178 * t199) * t139 - (t235 * t136 - t143 * t159 + t144 * t165) * t266) * t138 (t156 * t139 - (t235 * t135 - t143 * t153 + t144 * t164) * t266) * t138, 0, 0; ((-t153 * t221 - t167 * t226) * t147 - (-t153 * t226 + t167 * t221) * t264) * t142 ((t158 * t221 - t169 * t226) * t147 - (t158 * t226 + t169 * t221) * t264) * t142 ((t160 * t221 + t226 * t259) * t147 - (t160 * t226 - t221 * t259) * t264) * t142 (-t221 * t147 + t226 * t264) * t155 * t142, t249 * t142, 0;];
Ja_rot  = t1;
