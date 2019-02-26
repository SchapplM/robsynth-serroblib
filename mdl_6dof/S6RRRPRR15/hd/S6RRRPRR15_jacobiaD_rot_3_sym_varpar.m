% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRPRR15
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d5,d6]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRRPRR15_jacobiaD_rot_3_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR15_jacobiaD_rot_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR15_jacobiaD_rot_3_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:24:15
% EndTime: 2019-02-26 22:24:16
% DurationCPUTime: 0.83s
% Computational Cost: add. (2555->103), mult. (7918->233), div. (442->12), fcn. (10062->13), ass. (0->106)
t218 = sin(pkin(6));
t219 = cos(pkin(7));
t220 = cos(pkin(6));
t217 = sin(pkin(7));
t225 = cos(qJ(2));
t268 = t217 * t225;
t206 = -t218 * t268 + t219 * t220;
t203 = 0.1e1 / t206;
t223 = sin(qJ(1));
t260 = t223 * t225;
t222 = sin(qJ(2));
t226 = cos(qJ(1));
t262 = t222 * t226;
t238 = t220 * t262 + t260;
t267 = t218 * t222;
t204 = 0.1e1 / t206 ^ 2;
t259 = t225 * t226;
t261 = t223 * t222;
t207 = -t220 * t259 + t261;
t265 = t218 * t226;
t242 = -t207 * t217 + t219 * t265;
t272 = t242 * t204;
t284 = t217 * (t203 * t238 + t267 * t272);
t189 = atan2(t242, t206);
t184 = sin(t189);
t185 = cos(t189);
t170 = t184 * t242 + t185 * t206;
t167 = 0.1e1 / t170;
t221 = sin(qJ(3));
t224 = cos(qJ(3));
t237 = t220 * t261 - t259;
t239 = t220 * t260 + t262;
t266 = t218 * t223;
t252 = t217 * t266;
t240 = -t219 * t239 + t252;
t181 = t240 * t221 - t224 * t237;
t175 = 0.1e1 / t181;
t168 = 0.1e1 / t170 ^ 2;
t176 = 0.1e1 / t181 ^ 2;
t200 = -t217 * t239 - t219 * t266;
t197 = t200 ^ 2;
t163 = t168 * t197 + 0.1e1;
t192 = t207 * qJD(1) + t237 * qJD(2);
t258 = qJD(1) * t218;
t249 = t226 * t258;
t182 = t192 * t217 - t219 * t249;
t276 = t182 * t168;
t196 = t242 ^ 2;
t188 = t196 * t204 + 0.1e1;
t186 = 0.1e1 / t188;
t194 = t239 * qJD(1) + t238 * qJD(2);
t250 = t223 * t258;
t183 = -t194 * t217 - t219 * t250;
t257 = qJD(2) * t218;
t269 = t217 * t222;
t245 = t257 * t269;
t244 = t204 * t245;
t233 = t183 * t203 - t242 * t244;
t159 = t233 * t186;
t243 = -t184 * t206 + t185 * t242;
t155 = t243 * t159 + t184 * t183 + t185 * t245;
t282 = t155 * t167 * t168;
t283 = (-t197 * t282 + t200 * t276) / t163 ^ 2;
t193 = t238 * qJD(1) + t239 * qJD(2);
t235 = t192 * t219 + t217 * t249;
t165 = t181 * qJD(3) - t193 * t221 - t235 * t224;
t263 = t219 * t224;
t270 = t237 * t221;
t180 = -t224 * t252 + t239 * t263 - t270;
t174 = t180 ^ 2;
t173 = t174 * t176 + 0.1e1;
t277 = t176 * t180;
t166 = -t193 * t224 + t235 * t221 + (t240 * t224 + t270) * qJD(3);
t279 = t166 * t175 * t176;
t281 = (t165 * t277 - t174 * t279) / t173 ^ 2;
t205 = t203 * t204;
t280 = (-t196 * t205 * t245 + t183 * t272) / t188 ^ 2;
t278 = t168 * t200;
t275 = t184 * t200;
t274 = t185 * t200;
t273 = t242 * t203;
t271 = t238 * t221;
t264 = t219 * t221;
t256 = -0.2e1 * t283;
t255 = -0.2e1 * t282;
t254 = 0.2e1 * t281;
t253 = 0.2e1 * t280;
t251 = t217 * t265;
t248 = -0.2e1 * t203 * t280;
t247 = 0.2e1 * t180 * t279;
t246 = t217 * t250;
t241 = t207 * t219 + t251;
t190 = -t221 * t239 - t237 * t263;
t191 = -t224 * t239 + t237 * t264;
t234 = t184 + (t185 * t273 - t184) * t186;
t179 = t241 * t221 - t224 * t238;
t216 = t217 ^ 2;
t195 = t237 * qJD(1) + t207 * qJD(2);
t178 = -t241 * t224 - t271;
t171 = 0.1e1 / t173;
t161 = 0.1e1 / t163;
t160 = t186 * t284;
t158 = t234 * t200;
t156 = (-t184 * t238 + t185 * t267) * t217 - t243 * t160;
t154 = t253 * t284 + (t195 * t203 * t217 + (-t183 * t204 * t269 + (t204 * t238 * t216 * t222 + (0.2e1 * t205 * t216 * t218 * t222 ^ 2 - t204 * t268) * t242) * qJD(2)) * t218) * t186;
t1 = [t200 * t248 + (t182 * t203 - t200 * t244) * t186, t154, 0, 0, 0, 0; t242 * t167 * t256 + (t183 * t167 + (-t155 * t242 + t158 * t182) * t168) * t161 + ((t158 * t255 + t234 * t276) * t161 + (t158 * t256 + ((-t159 * t186 * t273 + t253) * t275 + (t242 * t248 + t159 + (-t159 + t233) * t186) * t274) * t161) * t168) * t200, 0.2e1 * (t167 * t217 * t237 - t156 * t278) * t283 + ((t243 * t154 - (-t170 * t159 + t183 * t185) * t160) * t278 + (t200 * t255 + t276) * t156 + (-t193 * t167 + (t237 * t155 + (-t159 * t238 + t225 * t257) * t274 + (t195 + (qJD(2) * t160 - t159) * t267) * t275) * t168) * t217) * t161, 0, 0, 0, 0; (-t175 * t178 + t179 * t277) * t254 + ((-t194 * t263 + t195 * t221 + t224 * t246) * t175 + t179 * t247 + (-t178 * t166 - (t194 * t264 + t195 * t224 - t221 * t246) * t180 - t179 * t165) * t176 + (t179 * t175 - (t207 * t263 + t224 * t251 + t271) * t277) * qJD(3)) * t171 (-t175 * t190 + t191 * t277) * t254 + ((t191 * qJD(3) + t192 * t221 - t193 * t263) * t175 + t191 * t247 + (-t190 * t166 - (-t190 * qJD(3) + t192 * t224 + t193 * t264) * t180 - t191 * t165) * t176) * t171, -0.2e1 * t281 + 0.2e1 * (t165 * t176 * t171 + (-t171 * t279 - t176 * t281) * t180) * t180, 0, 0, 0;];
JaD_rot  = t1;
