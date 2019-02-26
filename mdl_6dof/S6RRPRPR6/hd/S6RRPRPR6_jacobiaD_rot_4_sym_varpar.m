% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:40
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR6_jacobiaD_rot_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:40:37
% EndTime: 2019-02-26 21:40:38
% DurationCPUTime: 0.94s
% Computational Cost: add. (3306->95), mult. (9869->202), div. (448->12), fcn. (12731->13), ass. (0->96)
t216 = sin(qJ(4));
t218 = cos(qJ(4));
t271 = sin(pkin(11));
t273 = cos(pkin(6));
t241 = t273 * t271;
t272 = cos(pkin(11));
t242 = t273 * t272;
t274 = sin(qJ(2));
t275 = cos(qJ(2));
t206 = t275 * t241 + t274 * t242;
t208 = t274 * t271 - t275 * t272;
t217 = sin(qJ(1));
t219 = cos(qJ(1));
t238 = t217 * t206 + t219 * t208;
t215 = sin(pkin(6));
t256 = t215 * t217;
t233 = t216 * t238 + t218 * t256;
t278 = t233 * qJD(4);
t230 = t275 * t271 + t274 * t272;
t205 = t230 * t215;
t196 = qJD(2) * t205;
t204 = t208 * t215;
t201 = 0.1e1 / t204 ^ 2;
t258 = t196 * t201;
t277 = t274 * t241 - t275 * t242;
t229 = t208 * qJD(2);
t188 = -t217 * t230 - t219 * t277;
t176 = atan2(t188, t204);
t171 = sin(t176);
t172 = cos(t176);
t185 = t188 ^ 2;
t175 = t185 * t201 + 0.1e1;
t173 = 0.1e1 / t175;
t200 = 0.1e1 / t204;
t260 = t188 * t200;
t276 = (t172 * t260 - t171) * t173 + t171;
t160 = t171 * t188 + t172 * t204;
t157 = 0.1e1 / t160;
t184 = t216 * t256 - t218 * t238;
t178 = 0.1e1 / t184;
t158 = 0.1e1 / t160 ^ 2;
t179 = 0.1e1 / t184 ^ 2;
t226 = t217 * t277;
t191 = -t219 * t230 + t226;
t186 = t191 ^ 2;
t156 = t186 * t158 + 0.1e1;
t199 = t206 * qJD(2);
t166 = t188 * qJD(1) - t217 * t199 - t219 * t229;
t264 = t166 * t158;
t254 = qJD(1) * t219;
t169 = qJD(1) * t226 - t219 * t199 + t217 * t229 - t230 * t254;
t237 = t169 * t200 - t188 * t258;
t151 = t237 * t173;
t240 = -t171 * t204 + t172 * t188;
t147 = t240 * t151 + t171 * t169 + t172 * t196;
t269 = t147 * t157 * t158;
t270 = (-t186 * t269 - t191 * t264) / t156 ^ 2;
t239 = -t219 * t206 + t217 * t208;
t259 = t188 * t205;
t235 = -t200 * t239 + t201 * t259;
t152 = t235 * t173;
t148 = -t240 * t152 + t171 * t239 + t172 * t205;
t268 = t148 * t191;
t198 = t277 * qJD(2);
t207 = t230 * qJD(2);
t167 = t239 * qJD(1) + t217 * t198 - t219 * t207;
t248 = t215 * t254;
t161 = t184 * qJD(4) + t167 * t216 - t218 * t248;
t177 = t233 ^ 2;
t165 = t177 * t179 + 0.1e1;
t261 = t179 * t233;
t162 = t167 * t218 + t216 * t248 + t278;
t265 = t162 * t178 * t179;
t267 = (-t161 * t261 - t177 * t265) / t165 ^ 2;
t257 = t200 * t258;
t266 = (t188 * t201 * t169 - t185 * t257) / t175 ^ 2;
t263 = t171 * t191;
t262 = t172 * t191;
t255 = t215 * t219;
t253 = -0.2e1 * t270;
t252 = -0.2e1 * t269;
t251 = 0.2e1 * t267;
t250 = 0.2e1 * t266;
t249 = qJD(1) * t256;
t247 = -0.2e1 * t200 * t266;
t246 = -0.2e1 * t233 * t265;
t236 = -t216 * t178 - t218 * t261;
t234 = -t216 * t239 + t218 * t255;
t182 = t216 * t255 + t218 * t239;
t228 = t238 * qJD(1) + t219 * t198 + t217 * t207;
t197 = t215 * t229;
t163 = 0.1e1 / t165;
t154 = 0.1e1 / t156;
t150 = t276 * t191;
t146 = t235 * t250 + (0.2e1 * t257 * t259 + t228 * t200 + (-t169 * t205 + t188 * t197 - t196 * t239) * t201) * t173;
t1 = [t191 * t247 + (-t166 * t200 - t191 * t258) * t173, t146, 0, 0, 0, 0; t188 * t157 * t253 + (t169 * t157 + (-t147 * t188 - t150 * t166) * t158) * t154 + ((t150 * t252 - t276 * t264) * t154 + (t150 * t253 + ((-t151 * t173 * t260 + t250) * t263 + (t188 * t247 + t151 + (-t151 + t237) * t173) * t262) * t154) * t158) * t191, 0.2e1 * (t157 * t238 - t158 * t268) * t270 + (t167 * t157 + t252 * t268 + (t238 * t147 - t148 * t166 + (t146 * t188 - t152 * t169 - t197 + (t152 * t204 + t239) * t151) * t262 + (-t146 * t204 + t152 * t196 + t228 + (t152 * t188 - t205) * t151) * t263) * t158) * t154, 0, 0, 0, 0; (t178 * t234 - t182 * t261) * t251 + ((t182 * qJD(4) + t216 * t228 + t218 * t249) * t178 + t182 * t246 + (t234 * t162 + (t234 * qJD(4) - t216 * t249 + t218 * t228) * t233 - t182 * t161) * t179) * t163, t236 * t191 * t251 + (t236 * t166 + ((qJD(4) * t178 + t246) * t218 + (-t161 * t218 + (-t162 - t278) * t216) * t179) * t191) * t163, 0, -0.2e1 * t267 - 0.2e1 * (t161 * t179 * t163 - (-t163 * t265 - t179 * t267) * t233) * t233, 0, 0;];
JaD_rot  = t1;
