% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6PPRRRR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_jacobiaD_rot_5_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:42:43
% EndTime: 2019-02-26 19:42:44
% DurationCPUTime: 0.67s
% Computational Cost: add. (3188->68), mult. (9016->146), div. (299->12), fcn. (11669->15), ass. (0->82)
t242 = sin(pkin(7));
t243 = sin(pkin(6));
t245 = cos(pkin(6));
t241 = sin(pkin(13));
t283 = sin(pkin(12));
t267 = t283 * t241;
t244 = cos(pkin(13));
t284 = cos(pkin(12));
t268 = t284 * t244;
t285 = cos(pkin(7));
t289 = (t245 * t268 - t267) * t285 - t242 * t243 * t284;
t266 = t283 * t244;
t269 = t284 * t241;
t256 = t245 * t266 + t269;
t271 = t243 * t283;
t288 = -t242 * t271 + t256 * t285;
t232 = t245 * t269 + t266;
t246 = sin(qJ(3));
t286 = cos(qJ(3));
t218 = t232 * t286 + t289 * t246;
t270 = t244 * t285;
t273 = t242 * t245;
t287 = (-t241 * t246 + t286 * t270) * t243 + t286 * t273;
t216 = t232 * t246 - t289 * t286;
t209 = atan2(-t216, -t287);
t204 = sin(t209);
t205 = cos(t209);
t192 = -t204 * t216 - t205 * t287;
t189 = 0.1e1 / t192;
t233 = -t245 * t267 + t268;
t220 = t233 * t286 - t288 * t246;
t229 = t256 * t242 + t285 * t271;
t240 = qJ(4) + qJ(5);
t237 = sin(t240);
t238 = cos(t240);
t203 = t220 * t238 + t229 * t237;
t199 = 0.1e1 / t203;
t224 = 0.1e1 / t287;
t190 = 0.1e1 / t192 ^ 2;
t200 = 0.1e1 / t203 ^ 2;
t225 = 0.1e1 / t287 ^ 2;
t214 = t216 ^ 2;
t208 = t214 * t225 + 0.1e1;
t206 = 0.1e1 / t208;
t211 = t218 * qJD(3);
t228 = t246 * t273 + (t286 * t241 + t246 * t270) * t243;
t223 = t228 * qJD(3);
t277 = t216 * t225;
t183 = (t211 * t224 + t223 * t277) * t206;
t261 = t204 * t287 - t205 * t216;
t180 = t261 * t183 - t204 * t211 + t205 * t223;
t282 = t180 * t189 * t190;
t202 = t220 * t237 - t229 * t238;
t198 = t202 ^ 2;
t195 = t198 * t200 + 0.1e1;
t219 = t233 * t246 + t288 * t286;
t212 = t219 * qJD(3);
t239 = qJD(4) + qJD(5);
t265 = t229 * t239 - t212;
t275 = t220 * t239;
t196 = t265 * t237 + t238 * t275;
t278 = t200 * t202;
t197 = -t237 * t275 + t265 * t238;
t279 = t197 * t199 * t200;
t281 = (t196 * t278 - t198 * t279) / t195 ^ 2;
t280 = t190 * t219;
t276 = t216 * t228;
t274 = t223 * t224 * t225;
t272 = -0.2e1 * t281;
t259 = -t199 * t237 + t238 * t278;
t258 = t218 * t224 + t225 * t276;
t222 = t287 * qJD(3);
t215 = t219 ^ 2;
t213 = t220 * qJD(3);
t210 = t216 * qJD(3);
t193 = 0.1e1 / t195;
t187 = t215 * t190 + 0.1e1;
t184 = t258 * t206;
t181 = t261 * t184 - t204 * t218 + t205 * t228;
t179 = -0.2e1 * t258 / t208 ^ 2 * (t211 * t277 + t214 * t274) + (0.2e1 * t274 * t276 - t210 * t224 + (t211 * t228 + t216 * t222 + t218 * t223) * t225) * t206;
t177 = t272 + 0.2e1 * (t193 * t196 * t200 + (-t193 * t279 - t200 * t281) * t202) * t202;
t1 = [0, 0, t179, 0, 0, 0; 0, 0, 0.2e1 * (t181 * t280 - t189 * t220) / t187 ^ 2 * (t213 * t280 - t215 * t282) + (-t212 * t189 + (-t220 * t180 - t181 * t213) * t190 + (0.2e1 * t181 * t282 + (-(-t179 * t216 - t184 * t211 + t222 + (t184 * t287 - t218) * t183) * t205 - (t179 * t287 - t184 * t223 + t210 + (t184 * t216 - t228) * t183) * t204) * t190) * t219) / t187, 0, 0, 0; 0, 0, t259 * t219 * t272 + (t259 * t213 + ((-t199 * t239 - 0.2e1 * t202 * t279) * t238 + (t196 * t238 + (-t202 * t239 + t197) * t237) * t200) * t219) * t193, t177, t177, 0;];
JaD_rot  = t1;
