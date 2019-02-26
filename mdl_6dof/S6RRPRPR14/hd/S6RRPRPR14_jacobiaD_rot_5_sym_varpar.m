% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR14
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR14_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR14_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR14_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:26
% EndTime: 2019-02-26 21:45:27
% DurationCPUTime: 1.26s
% Computational Cost: add. (3617->117), mult. (10934->241), div. (683->12), fcn. (13943->11), ass. (0->106)
t196 = cos(qJ(4));
t192 = cos(pkin(6));
t197 = cos(qJ(2));
t198 = cos(qJ(1));
t235 = t197 * t198;
t194 = sin(qJ(2));
t195 = sin(qJ(1));
t238 = t194 * t195;
t213 = t192 * t235 - t238;
t193 = sin(qJ(4));
t191 = sin(pkin(6));
t239 = t191 * t198;
t225 = t193 * t239;
t214 = -t196 * t213 + t225;
t167 = t214 ^ 2;
t240 = t191 * t197;
t184 = t192 * t193 + t196 * t240;
t179 = 0.1e1 / t184 ^ 2;
t158 = t167 * t179 + 0.1e1;
t154 = 0.1e1 / t158;
t222 = t196 * t239;
t189 = qJD(4) * t222;
t236 = t195 * t197;
t237 = t194 * t198;
t212 = t192 * t236 + t237;
t242 = t191 * t195;
t187 = t192 * t237 + t236;
t264 = t187 * qJD(2);
t148 = -(qJD(1) * t212 + t264) * t196 - t189 + (qJD(1) * t242 - qJD(4) * t213) * t193;
t185 = t192 * t196 - t193 * t240;
t243 = t191 * t194;
t221 = qJD(2) * t243;
t169 = qJD(4) * t185 - t196 * t221;
t178 = 0.1e1 / t184;
t246 = t214 * t179;
t217 = -t148 * t178 - t169 * t246;
t136 = t217 * t154;
t159 = atan2(t214, t184);
t152 = sin(t159);
t153 = cos(t159);
t218 = -t152 * t184 + t153 * t214;
t132 = t136 * t218 - t148 * t152 + t153 * t169;
t147 = t152 * t214 + t153 * t184;
t145 = 0.1e1 / t147 ^ 2;
t271 = t132 * t145;
t144 = 0.1e1 / t147;
t270 = t144 * t271;
t241 = t191 * t196;
t171 = t193 * t212 + t195 * t241;
t211 = t192 * t238 - t235;
t205 = t211 * qJD(2);
t268 = -qJD(1) * t213 + t205;
t150 = qJD(1) * t225 + t171 * qJD(4) + t196 * t268;
t261 = -t193 * t242 + t196 * t212;
t260 = -0.2e1 * t261;
t220 = t260 * t270;
t269 = -t145 * t150 + t220;
t215 = t193 * t213 + t222;
t149 = qJD(1) * t171 + t214 * qJD(4) + t193 * t264;
t163 = -qJD(1) * t187 - qJD(2) * t212;
t182 = 0.1e1 / t211 ^ 2;
t251 = t163 * t182;
t267 = t169 * t179;
t226 = t214 * t243;
t210 = t178 * t187 + t179 * t226;
t266 = t196 * t210;
t181 = 0.1e1 / t211;
t165 = t261 ^ 2;
t143 = t145 * t165 + 0.1e1;
t254 = t145 * t261;
t259 = (-t150 * t254 - t165 * t270) / t143 ^ 2;
t151 = qJD(1) * t215 + qJD(4) * t261 - t193 * t205;
t166 = t171 ^ 2;
t160 = t166 * t182 + 0.1e1;
t248 = t171 * t182;
t250 = t181 * t251;
t257 = (t151 * t248 + t166 * t250) / t160 ^ 2;
t249 = t178 * t267;
t256 = (-t148 * t246 - t167 * t249) / t158 ^ 2;
t253 = t152 * t261;
t252 = t153 * t261;
t247 = t214 * t178;
t245 = t214 * t185;
t244 = t211 * t196;
t234 = qJD(2) * t197;
t233 = qJD(4) * t193;
t232 = 0.2e1 * t259;
t231 = -0.2e1 * t256;
t230 = 0.2e1 * t256;
t228 = -0.2e1 * t171 * t187;
t227 = t181 * t257;
t219 = t178 * t230;
t216 = -t178 * t215 + t179 * t245;
t209 = t182 * t212;
t208 = -t152 + (-t153 * t247 + t152) * t154;
t168 = -qJD(4) * t184 + t193 * t221;
t164 = -qJD(1) * t211 + qJD(2) * t213;
t156 = 0.1e1 / t160;
t141 = 0.1e1 / t143;
t140 = t154 * t266;
t139 = t216 * t154;
t134 = (t152 * t187 - t153 * t243) * t196 + t218 * t140;
t133 = -t139 * t218 + t152 * t215 + t153 * t185;
t131 = t216 * t230 + (0.2e1 * t245 * t249 - t149 * t178 + (t148 * t185 - t168 * t214 - t169 * t215) * t179) * t154;
t129 = t231 * t266 + (-t210 * t233 + (-0.2e1 * t226 * t249 + t164 * t178 + (-t169 * t187 + (-t148 * t194 + t214 * t234) * t191) * t179) * t196) * t154;
t1 = [-t261 * t219 + (-t150 * t178 - t261 * t267) * t154, t129, 0, t131, 0, 0; -0.2e1 * t214 * t144 * t259 + ((qJD(1) * t261 + t196 * t264 + t213 * t233 + t189) * t144 - t214 * t271 + (t208 * t150 - ((t136 * t154 * t247 + t231) * t152 + (t214 * t219 - t136 + (t136 - t217) * t154) * t153) * t261) * t254) * t141 - (t141 * t269 - t254 * t232) * t208 * t261 (-t134 * t254 - t144 * t244) * t232 + ((-t163 * t196 - t211 * t233) * t144 + t269 * t134 + (-t244 * t132 + (t129 * t214 - t140 * t148 + (t194 * t233 - t196 * t234) * t191 + (-t140 * t184 + t187 * t196) * t136) * t252 + (-t187 * t233 - t129 * t184 - t140 * t169 + t164 * t196 + (-t140 * t214 + t194 * t241) * t136) * t253) * t145) * t141, 0 (-t133 * t254 - t144 * t171) * t232 + (t133 * t220 + t151 * t144 + (-t171 * t132 - t133 * t150 + (t131 * t214 + t139 * t148 + t168 + (t139 * t184 + t215) * t136) * t252 + (-t131 * t184 + t139 * t169 - t149 + (t139 * t214 - t185) * t136) * t253) * t145) * t141, 0, 0; t182 * t228 * t257 + 0.2e1 * t215 * t227 + (t187 * t151 * t182 + t149 * t181 + t164 * t248 - t215 * t251 - t228 * t250) * t156, 0.2e1 * (-t171 * t209 - t193) * t257 + (qJD(4) * t196 + t151 * t209 + (-t182 * t268 + 0.2e1 * t212 * t250) * t171) * t156, 0, -t227 * t260 + (t150 * t181 - t251 * t261) * t156, 0, 0;];
JaD_rot  = t1;
