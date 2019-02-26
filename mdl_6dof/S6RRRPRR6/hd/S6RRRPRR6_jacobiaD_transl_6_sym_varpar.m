% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRR6
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:48
% EndTime: 2019-02-26 22:18:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (804->76), mult. (632->102), div. (0->0), fcn. (500->12), ass. (0->65)
t255 = qJ(3) + pkin(11);
t251 = qJ(5) + t255;
t245 = sin(t251);
t273 = -sin(qJ(3)) * pkin(3) - pkin(4) * sin(t255);
t254 = qJD(3) + qJD(5);
t300 = pkin(5) * t254;
t235 = t273 * qJD(3) - t245 * t300;
t246 = cos(t251);
t305 = cos(qJ(3)) * pkin(3) + pkin(4) * cos(t255);
t239 = pkin(5) * t246 + t305;
t237 = pkin(2) + t239;
t301 = pkin(5) * t245;
t238 = -t273 + t301;
t257 = sin(qJ(2));
t260 = cos(qJ(2));
t284 = t257 * qJD(4);
t295 = r_i_i_C(3) + pkin(10) + pkin(9) + qJ(4) + pkin(8);
t306 = t295 * t260;
t310 = (-t237 * t257 + t306) * qJD(2) + (pkin(7) + t238) * qJD(1) + t260 * t235 + t284;
t247 = qJ(6) + t251;
t242 = sin(t247);
t298 = r_i_i_C(2) * t242;
t243 = cos(t247);
t299 = r_i_i_C(1) * t243;
t267 = t237 - t298 + t299;
t263 = -t267 * t257 + t306;
t261 = cos(qJ(1));
t250 = qJD(6) + t254;
t276 = t250 * t260 - qJD(1);
t308 = t261 * t276;
t297 = r_i_i_C(2) * t243;
t272 = r_i_i_C(1) * t242 + t297;
t304 = t272 * t250 - t235;
t289 = qJD(1) * t260;
t275 = -t250 + t289;
t258 = sin(qJ(1));
t287 = qJD(2) * t257;
t281 = t258 * t287;
t303 = t275 * t261 - t281;
t293 = t250 * t257;
t285 = qJD(2) * t261;
t280 = t257 * t285;
t265 = t275 * t258 + t280;
t231 = t265 * t242 - t243 * t308;
t232 = t242 * t308 + t265 * t243;
t292 = t231 * r_i_i_C(1) + t232 * r_i_i_C(2);
t269 = t276 * t258;
t233 = t303 * t242 + t243 * t269;
t234 = t242 * t269 - t303 * t243;
t291 = -t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
t290 = qJD(1) * t258;
t288 = qJD(1) * t261;
t286 = qJD(2) * t260;
t283 = t246 * t300;
t282 = t250 * t299;
t279 = t295 * t257;
t274 = -t254 + t289;
t271 = t238 * t289 + t235;
t270 = t246 * (-t254 * t260 + qJD(1));
t236 = t305 * qJD(3) + t283;
t266 = qJD(1) * t239 - t236 * t260 + t238 * t287;
t264 = t236 + (-t237 * t260 - pkin(1) - t279) * qJD(1);
t262 = qJD(4) * t260 + t304 * t257 + (-t267 * t260 - t279) * qJD(2);
t240 = t293 * t298;
t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) - t310 * t258 + t264 * t261, t262 * t261 - t263 * t290, t271 * t258 + t266 * t261 + t292, -t257 * t290 + t260 * t285 (t261 * t270 + (t274 * t258 + t280) * t245) * pkin(5) + t292, t292; -t232 * r_i_i_C(1) + t231 * r_i_i_C(2) + t264 * t258 + t310 * t261, t262 * t258 + t263 * t288, t266 * t258 - t271 * t261 + t291, t257 * t288 + t258 * t286 (t258 * t270 + (-t261 * t274 + t281) * t245) * pkin(5) + t291, t291; 0, t263 * qJD(2) - t304 * t260 + t284, t240 + (-t236 - t282) * t257 + (-t238 - t272) * t286, t287, t240 + (-t282 - t283) * t257 + (-t272 - t301) * t286, -t286 * t297 + t240 + (-t242 * t286 - t243 * t293) * r_i_i_C(1);];
JaD_transl  = t1;
