% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:47:41
% EndTime: 2019-02-26 22:47:41
% DurationCPUTime: 0.30s
% Computational Cost: add. (576->71), mult. (524->104), div. (0->0), fcn. (386->10), ass. (0->64)
t272 = qJ(2) + qJ(3);
t269 = qJ(4) + t272;
t264 = sin(t269);
t276 = cos(qJ(5));
t331 = r_i_i_C(1) * t276 + pkin(4);
t292 = t331 * t264;
t273 = sin(qJ(5));
t313 = qJD(5) * t276;
t265 = cos(t269);
t270 = qJD(2) + qJD(3);
t266 = qJD(4) + t270;
t321 = t265 * t266;
t333 = t264 * t313 + t273 * t321;
t328 = pkin(10) + r_i_i_C(3);
t306 = t328 * t265;
t274 = sin(qJ(2));
t322 = pkin(2) * qJD(2);
t308 = t274 * t322;
t267 = sin(t272);
t326 = pkin(3) * t270;
t312 = t267 * t326;
t325 = pkin(4) * t264;
t332 = (t306 - t325) * t266 - t308 - t312;
t314 = qJD(5) * t273;
t299 = t264 * t314;
t329 = r_i_i_C(1) * t299 + t333 * r_i_i_C(2);
t268 = cos(t272);
t291 = t331 * t265;
t307 = t328 * t264;
t281 = (-t291 - t307) * t266 - t268 * t326;
t327 = pkin(3) * t267;
t323 = r_i_i_C(2) * t273;
t275 = sin(qJ(1));
t320 = t266 * t275;
t319 = t266 * t276;
t278 = cos(qJ(1));
t318 = t266 * t278;
t317 = t276 * t278;
t316 = qJD(1) * t275;
t315 = qJD(1) * t278;
t310 = t264 * t323;
t309 = qJD(1) * t323;
t305 = t328 * t275;
t304 = t264 * t319;
t294 = qJD(5) * t265 - qJD(1);
t293 = qJD(1) * t265 - qJD(5);
t290 = t331 * t278;
t289 = t329 * t278 + t316 * t292;
t288 = t294 * t273;
t287 = t278 * t264 * t309 + t329 * t275 + t315 * t306;
t286 = -t306 - t310;
t277 = cos(qJ(2));
t285 = -t277 * pkin(2) - pkin(3) * t268 - pkin(4) * t265 - pkin(1) - t307;
t284 = t264 * t318 + t293 * t275;
t282 = -t277 * t322 + t281;
t280 = -t265 * r_i_i_C(2) * t313 + (-t265 * t314 - t304) * r_i_i_C(1) + t328 * t321 + (-t325 + t310) * t266;
t279 = t280 - t312;
t271 = -pkin(9) - pkin(8) - pkin(7);
t263 = -t274 * pkin(2) - t327;
t245 = -t293 * t317 + (t288 + t304) * t275;
t244 = t294 * t276 * t275 + (-t264 * t320 + t293 * t278) * t273;
t243 = t284 * t276 + t278 * t288;
t242 = t284 * t273 - t294 * t317;
t1 = [t245 * r_i_i_C(1) + t244 * r_i_i_C(2) - t332 * t275 + (t271 * t275 + t285 * t278) * qJD(1) (-t263 + t286) * t316 + t282 * t278 + t289 (t286 + t327) * t316 + t281 * t278 + t289 (-t275 * t309 - t328 * t318) * t264 + (-qJD(1) * t305 - t266 * t290) * t265 + t289, t242 * r_i_i_C(1) + t243 * r_i_i_C(2), 0; -t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t332 * t278 + (-t271 * t278 + t285 * t275) * qJD(1) (t263 - t292) * t315 + t282 * t275 + t287 (-t292 - t327) * t315 + t281 * t275 + t287, -t291 * t320 + (-qJD(1) * t290 - t266 * t305) * t264 + t287, -t244 * r_i_i_C(1) + t245 * r_i_i_C(2), 0; 0, t279 - t308, t279, t280 (-t265 * t319 + t299) * r_i_i_C(2) - t333 * r_i_i_C(1), 0;];
JaD_transl  = t1;
