% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPPR10
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPPR10_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR10_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR10_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:57
% EndTime: 2019-02-26 22:08:58
% DurationCPUTime: 1.04s
% Computational Cost: add. (3640->122), mult. (11003->253), div. (691->12), fcn. (14028->11), ass. (0->105)
t197 = cos(pkin(6));
t200 = sin(qJ(1));
t202 = cos(qJ(2));
t258 = cos(qJ(1));
t223 = t258 * qJD(1);
t224 = qJD(2) * t258;
t196 = sin(pkin(6));
t229 = t196 * t258;
t199 = sin(qJ(2));
t242 = t200 * t199;
t230 = t197 * t242;
t240 = qJD(2) * t199;
t263 = -qJD(1) * t230 - t200 * t240 + (t197 * t224 + t223) * t202 - qJD(3) * t229;
t198 = sin(qJ(3));
t201 = cos(qJ(3));
t228 = t258 * t199;
t241 = t200 * t202;
t210 = -t197 * t228 - t241;
t167 = -t198 * t210 + t201 * t229;
t180 = t196 * t199 * t198 - t197 * t201;
t156 = atan2(-t167, t180);
t149 = sin(t156);
t150 = cos(t156);
t144 = -t149 * t167 + t150 * t180;
t142 = 0.1e1 / t144 ^ 2;
t227 = t258 * t202;
t185 = t227 - t230;
t244 = t196 * t201;
t214 = -t185 * t198 + t200 * t244;
t163 = t214 ^ 2;
t140 = t163 * t142 + 0.1e1;
t211 = -t197 * t241 - t228;
t159 = t210 * qJD(1) + t211 * qJD(2);
t245 = t196 * t200;
t173 = t185 * t201 + t198 * t245;
t219 = t196 * t223;
t145 = t173 * qJD(3) + t159 * t198 - t201 * t219;
t253 = t142 * t214;
t162 = t167 ^ 2;
t175 = 0.1e1 / t180 ^ 2;
t155 = t162 * t175 + 0.1e1;
t151 = 0.1e1 / t155;
t226 = qJD(1) * t245;
t239 = qJD(3) * t201;
t147 = t263 * t198 - t201 * t226 - t210 * t239;
t181 = t197 * t198 + t199 * t244;
t243 = t196 * t202;
t225 = qJD(2) * t243;
t165 = t181 * qJD(3) + t198 * t225;
t174 = 0.1e1 / t180;
t248 = t167 * t175;
t216 = -t147 * t174 + t165 * t248;
t133 = t216 * t151;
t217 = -t149 * t180 - t150 * t167;
t129 = t217 * t133 - t149 * t147 + t150 * t165;
t141 = 0.1e1 / t144;
t143 = t141 * t142;
t256 = t129 * t143;
t238 = 0.2e1 * (-t145 * t253 - t163 * t256) / t140 ^ 2;
t262 = t165 * t175;
t220 = t197 * t227;
t182 = t220 - t242;
t212 = -t174 * t182 + t243 * t248;
t261 = t198 * t212;
t148 = (qJD(3) * t210 + t226) * t198 + t263 * t201;
t177 = 0.1e1 / t211;
t178 = 0.1e1 / t211 ^ 2;
t260 = -0.2e1 * t167;
t259 = -0.2e1 * t214;
t250 = t174 * t262;
t255 = (t147 * t248 - t162 * t250) / t155 ^ 2;
t254 = t142 * t145;
t252 = t149 * t214;
t251 = t150 * t214;
t249 = t167 * t174;
t247 = t173 * t178;
t246 = t211 * t198;
t146 = t214 * qJD(3) + t159 * t201 + t198 * t219;
t164 = t173 ^ 2;
t157 = t164 * t178 + 0.1e1;
t158 = -qJD(1) * t220 - t202 * t224 + (qJD(2) * t197 + qJD(1)) * t242;
t179 = t177 * t178;
t237 = 0.2e1 * (-t164 * t179 * t158 + t146 * t247) / t157 ^ 2;
t236 = -0.2e1 * t255;
t235 = t143 * t259;
t234 = 0.2e1 * t173 * t179;
t233 = t174 * t255;
t232 = t142 * t252;
t231 = t142 * t251;
t222 = t250 * t260;
t169 = -t198 * t229 - t201 * t210;
t215 = -t169 * t174 + t181 * t248;
t209 = -t149 + (t150 * t249 + t149) * t151;
t166 = -t180 * qJD(3) + t201 * t225;
t160 = t211 * qJD(1) + t210 * qJD(2);
t153 = 0.1e1 / t157;
t138 = 0.1e1 / t140;
t137 = t151 * t261;
t136 = t215 * t151;
t132 = t209 * t214;
t131 = (-t149 * t182 + t150 * t243) * t198 + t217 * t137;
t130 = t217 * t136 - t149 * t169 + t150 * t181;
t128 = t215 * t236 + (t181 * t222 - t148 * t174 + (t147 * t181 + t165 * t169 + t166 * t167) * t175) * t151;
t126 = t236 * t261 + (t212 * t239 + (t222 * t243 - t160 * t174 + (t165 * t182 + (t147 * t202 - t167 * t240) * t196) * t175) * t198) * t151;
t1 = [t233 * t259 + (-t145 * t174 - t214 * t262) * t151, t126, t128, 0, 0, 0; t167 * t141 * t238 + (-t147 * t141 + (t129 * t167 + t132 * t145) * t142) * t138 - (-t132 * t142 * t238 + (-0.2e1 * t132 * t256 + (-t133 * t151 * t249 + t236) * t232 + (t233 * t260 - t133 + (t133 - t216) * t151) * t231 - t209 * t254) * t138) * t214 (-t131 * t253 - t141 * t246) * t238 + (-t131 * t254 + (t158 * t198 + t211 * t239) * t141 + (t131 * t235 - t142 * t246) * t129 + (-t126 * t167 - t137 * t147 + (-t198 * t240 + t202 * t239) * t196 + (-t137 * t180 - t182 * t198) * t133) * t231 + (-t182 * t239 - t126 * t180 - t137 * t165 - t160 * t198 + (t137 * t167 - t198 * t243) * t133) * t232) * t138 (-t130 * t253 - t141 * t173) * t238 + (t130 * t129 * t235 + t146 * t141 + (-t173 * t129 - t130 * t145 + (-t128 * t167 - t136 * t147 + t166 + (-t136 * t180 - t169) * t133) * t251 + (-t128 * t180 - t136 * t165 - t148 + (t136 * t167 - t181) * t133) * t252) * t142) * t138, 0, 0, 0; (-t169 * t177 + t182 * t247) * t237 + (t148 * t177 + t182 * t158 * t234 + (-t182 * t146 - t158 * t169 - t160 * t173) * t178) * t153 (t177 * t201 * t211 + t185 * t247) * t237 + (qJD(3) * t177 * t246 + (-t146 * t185 - t159 * t173) * t178 + (t185 * t234 + (t178 * t211 - t177) * t201) * t158) * t153, t214 * t177 * t237 + (t158 * t178 * t214 + t145 * t177) * t153, 0, 0, 0;];
JaD_rot  = t1;
