% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRRPR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:07
% EndTime: 2019-02-26 22:34:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (523->86), mult. (804->130), div. (0->0), fcn. (741->12), ass. (0->63)
t301 = sin(qJ(1));
t298 = cos(pkin(6));
t314 = qJD(2) * t298 + qJD(1);
t300 = sin(qJ(2));
t330 = t301 * t300;
t317 = t298 * t330;
t321 = qJD(2) * t300;
t303 = cos(qJ(2));
t304 = cos(qJ(1));
t327 = t304 * t303;
t271 = -qJD(1) * t317 - t301 * t321 + t314 * t327;
t295 = qJD(3) + qJD(4);
t297 = sin(pkin(6));
t331 = t297 * t304;
t311 = t295 * t331 - t271;
t296 = qJ(3) + qJ(4);
t291 = sin(t296);
t299 = sin(qJ(3));
t335 = pkin(3) * qJD(3);
t338 = pkin(4) * t295;
t272 = -t291 * t338 - t299 * t335;
t290 = pkin(12) + t296;
t287 = sin(t290);
t288 = cos(t290);
t313 = r_i_i_C(1) * t287 + r_i_i_C(2) * t288;
t339 = t313 * t295 - t272;
t280 = t299 * pkin(3) + pkin(4) * t291;
t337 = pkin(8) + t280;
t336 = r_i_i_C(3) + qJ(5) + pkin(10) + pkin(9);
t328 = t304 * t300;
t329 = t301 * t303;
t276 = t298 * t328 + t329;
t334 = t276 * t295;
t333 = t287 * t295;
t332 = t297 * t301;
t278 = -t317 + t327;
t322 = qJD(1) * t304;
t307 = -t278 * t295 + t297 * t322;
t309 = t298 * t329 + t328;
t269 = t276 * qJD(1) + t309 * qJD(2);
t312 = t295 * t332 - t269;
t264 = -t312 * t287 + t307 * t288;
t265 = t307 * t287 + t312 * t288;
t326 = t264 * r_i_i_C(1) - t265 * r_i_i_C(2);
t323 = qJD(1) * t301;
t308 = t297 * t323 - t334;
t315 = t311 * t288;
t325 = (t311 * t287 + t308 * t288) * r_i_i_C(1) + (-t308 * t287 + t315) * r_i_i_C(2);
t320 = qJD(2) * t303;
t306 = -t295 * t298 - t297 * t320;
t319 = t295 * t297 * t300;
t324 = (t306 * t287 - t288 * t319) * r_i_i_C(1) + (t287 * t319 + t306 * t288) * r_i_i_C(2);
t292 = cos(t296);
t302 = cos(qJ(3));
t281 = t302 * pkin(3) + pkin(4) * t292;
t316 = t298 * t327;
t279 = pkin(2) + t281;
t310 = t288 * r_i_i_C(1) - t287 * r_i_i_C(2) + t279;
t275 = -t316 + t330;
t273 = t292 * t338 + t302 * t335;
t270 = t309 * qJD(1) + t276 * qJD(2);
t268 = -qJD(1) * t316 - t304 * t320 + t314 * t330;
t1 = [(t276 * t333 + t315) * r_i_i_C(1) + (t271 * t287 + t288 * t334) * r_i_i_C(2) - t271 * t279 - t276 * t272 - t275 * qJD(5) - pkin(1) * t322 - t336 * t270 + ((-r_i_i_C(2) * t333 + t273) * t304 + (-t313 - t337) * t323) * t297, t278 * qJD(5) + t310 * t268 - t336 * t269 + t309 * t339, t269 * t280 - t278 * t273 + (t272 * t301 + t281 * t322) * t297 + t326 (-t312 * t291 + t307 * t292) * pkin(4) + t326, -t268, 0; t273 * t332 + t265 * r_i_i_C(1) + t264 * r_i_i_C(2) + t309 * qJD(5) - t269 * t279 + t278 * t272 - t336 * t268 + (-pkin(1) * t301 + t337 * t331) * qJD(1), t276 * qJD(5) - t310 * t270 + t336 * t271 + t339 * t275, -t271 * t280 - t276 * t273 + (-t272 * t304 + t281 * t323) * t297 + t325 (t311 * t291 + t308 * t292) * pkin(4) + t325, t270, 0; 0 (qJD(5) * t300 - t339 * t303 + (-t310 * t300 + t336 * t303) * qJD(2)) * t297, t298 * t272 + (-t273 * t300 - t280 * t320) * t297 + t324 (t306 * t291 - t292 * t319) * pkin(4) + t324, t297 * t321, 0;];
JaD_transl  = t1;
