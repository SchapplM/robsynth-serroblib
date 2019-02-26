% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PRRRPR2_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPR2_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRRPR2_jacobiaD_rot_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:11:07
% EndTime: 2019-02-26 20:11:08
% DurationCPUTime: 1.07s
% Computational Cost: add. (8455->106), mult. (12233->219), div. (792->12), fcn. (15777->13), ass. (0->109)
t235 = sin(pkin(11));
t238 = cos(pkin(11));
t240 = sin(qJ(2));
t239 = cos(pkin(6));
t241 = cos(qJ(2));
t265 = t239 * t241;
t222 = -t235 * t240 + t238 * t265;
t218 = t222 * qJD(2);
t266 = t239 * t240;
t223 = t235 * t241 + t238 * t266;
t233 = qJ(3) + qJ(4);
t230 = sin(t233);
t232 = qJD(3) + qJD(4);
t236 = sin(pkin(6));
t269 = t236 * t238;
t255 = t230 * t269;
t231 = cos(t233);
t271 = t231 * t232;
t186 = t218 * t230 + t223 * t271 - t232 * t255;
t208 = t223 * t230 + t231 * t269;
t206 = t208 ^ 2;
t268 = t236 * t240;
t257 = t230 * t268;
t216 = -t239 * t231 + t257;
t214 = 0.1e1 / t216 ^ 2;
t200 = t206 * t214 + 0.1e1;
t198 = 0.1e1 / t200;
t263 = qJD(2) * t241;
t249 = t232 * t239 + t236 * t263;
t256 = t231 * t268;
t204 = t230 * t249 + t232 * t256;
t213 = 0.1e1 / t216;
t276 = t208 * t214;
t170 = (-t186 * t213 + t204 * t276) * t198;
t201 = atan2(-t208, t216);
t196 = sin(t201);
t197 = cos(t201);
t252 = -t196 * t216 - t197 * t208;
t166 = t170 * t252 - t196 * t186 + t197 * t204;
t180 = -t196 * t208 + t197 * t216;
t177 = 0.1e1 / t180;
t178 = 0.1e1 / t180 ^ 2;
t289 = t166 * t177 * t178;
t267 = t236 * t241;
t248 = -t213 * t222 + t267 * t276;
t288 = t230 * t248;
t258 = t235 * t266;
t225 = t238 * t241 - t258;
t270 = t235 * t236;
t211 = t225 * t230 - t231 * t270;
t287 = 0.2e1 * t211 * t289;
t277 = t204 * t213 * t214;
t286 = -0.2e1 * (t186 * t276 - t206 * t277) / t200 ^ 2;
t212 = t225 * t231 + t230 * t270;
t237 = cos(pkin(12));
t224 = t235 * t265 + t238 * t240;
t234 = sin(pkin(12));
t274 = t224 * t234;
t195 = t212 * t237 + t274;
t191 = 0.1e1 / t195;
t192 = 0.1e1 / t195 ^ 2;
t285 = t178 * t211;
t220 = t224 * qJD(2);
t253 = t232 * t270 - t220;
t272 = t230 * t232;
t189 = -t225 * t272 + t231 * t253;
t221 = -qJD(2) * t258 + t238 * t263;
t185 = t189 * t237 + t221 * t234;
t284 = t185 * t191 * t192;
t188 = t225 * t271 + t230 * t253;
t283 = t188 * t178;
t282 = t191 * t234;
t273 = t224 * t237;
t194 = t212 * t234 - t273;
t281 = t192 * t194;
t280 = t194 * t237;
t279 = t196 * t211;
t278 = t197 * t211;
t275 = t224 * t230;
t264 = qJD(2) * t240;
t207 = t211 ^ 2;
t176 = t207 * t178 + 0.1e1;
t262 = 0.2e1 * (-t207 * t289 + t211 * t283) / t176 ^ 2;
t190 = t194 ^ 2;
t183 = t190 * t192 + 0.1e1;
t184 = t189 * t234 - t221 * t237;
t261 = 0.2e1 * (t184 * t281 - t190 * t284) / t183 ^ 2;
t259 = t194 * t284;
t254 = -0.2e1 * t208 * t277;
t210 = t223 * t231 - t255;
t217 = t239 * t230 + t256;
t251 = -t210 * t213 + t217 * t276;
t250 = -t221 * t231 + t224 * t272;
t219 = t223 * qJD(2);
t205 = t231 * t249 - t232 * t257;
t203 = t225 * t234 - t231 * t273;
t202 = -t225 * t237 - t231 * t274;
t187 = -t223 * t272 + (-t232 * t269 + t218) * t231;
t181 = 0.1e1 / t183;
t174 = 0.1e1 / t176;
t172 = t198 * t288;
t171 = t251 * t198;
t168 = (-t196 * t222 + t197 * t267) * t230 + t252 * t172;
t167 = t171 * t252 - t196 * t210 + t197 * t217;
t165 = t251 * t286 + (t217 * t254 - t187 * t213 + (t186 * t217 + t204 * t210 + t205 * t208) * t214) * t198;
t163 = t286 * t288 + (t248 * t271 + (t254 * t267 + t213 * t219 + (t204 * t222 + (t186 * t241 - t208 * t264) * t236) * t214) * t230) * t198;
t162 = (-t192 * t280 + t282) * t211 * t261 + (-0.2e1 * t211 * t237 * t259 - t188 * t282 + (t188 * t280 + (t184 * t237 + t185 * t234) * t211) * t192) * t181;
t161 = (t167 * t285 - t177 * t212) * t262 + (t167 * t287 + t189 * t177 + (-t212 * t166 - t167 * t188 - (-t165 * t208 - t171 * t186 + t205 + (-t171 * t216 - t210) * t170) * t278 - (-t165 * t216 - t171 * t204 - t187 + (t171 * t208 - t217) * t170) * t279) * t178) * t174;
t1 = [0, t163, t165, t165, 0, 0; 0 (t168 * t285 + t177 * t275) * t262 + ((-t221 * t230 - t224 * t271) * t177 + (-t283 + t287) * t168 + (t275 * t166 - (-t163 * t208 - t172 * t186 + (-t230 * t264 + t241 * t271) * t236 + (-t172 * t216 - t222 * t230) * t170) * t278 - (-t222 * t271 - t163 * t216 - t172 * t204 + t219 * t230 + (t172 * t208 - t230 * t267) * t170) * t279) * t178) * t174, t161, t161, 0, 0; 0 (-t191 * t202 + t203 * t281) * t261 + ((t220 * t237 + t234 * t250) * t191 + 0.2e1 * t203 * t259 + (-t202 * t185 - (-t220 * t234 + t237 * t250) * t194 - t203 * t184) * t192) * t181, t162, t162, 0, 0;];
JaD_rot  = t1;
