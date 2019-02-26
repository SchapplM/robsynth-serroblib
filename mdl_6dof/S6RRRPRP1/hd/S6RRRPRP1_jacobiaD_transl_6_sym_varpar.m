% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPRP1
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:09
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP1_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:09:31
% EndTime: 2019-02-26 22:09:31
% DurationCPUTime: 0.38s
% Computational Cost: add. (538->79), mult. (547->106), div. (0->0), fcn. (416->10), ass. (0->64)
t338 = pkin(5) + r_i_i_C(1);
t275 = qJD(2) + qJD(3);
t278 = sin(qJ(5));
t329 = r_i_i_C(2) * t278;
t345 = t275 * t329 + qJD(6);
t276 = qJ(2) + qJ(3);
t271 = pkin(10) + t276;
t268 = sin(t271);
t269 = cos(t271);
t281 = cos(qJ(5));
t317 = qJD(5) * t281;
t318 = qJD(5) * t278;
t344 = t345 * t269 + (r_i_i_C(2) * t317 + t318 * t338) * t268;
t280 = sin(qJ(1));
t283 = cos(qJ(1));
t301 = qJD(1) * t269 - qJD(5);
t295 = t301 * t283;
t302 = qJD(5) * t269 - qJD(1);
t297 = t302 * t281;
t322 = t275 * t280;
t312 = t268 * t322;
t242 = t280 * t297 + (t295 - t312) * t278;
t270 = pkin(5) * t281 + pkin(4);
t343 = r_i_i_C(1) * t281 + t270;
t272 = sin(t276);
t279 = sin(qJ(2));
t326 = pkin(2) * qJD(2);
t313 = t279 * t326;
t277 = -qJ(6) - pkin(9);
t327 = r_i_i_C(3) - t277;
t331 = pkin(3) * t275;
t339 = (-pkin(5) * t318 + t327 * t275) * t269 + (pkin(5) * t278 + pkin(7) + pkin(8) + qJ(4)) * qJD(1) - (t270 * t275 - qJD(6)) * t268 - t272 * t331 - t313;
t333 = pkin(3) * t272;
t273 = cos(t276);
t332 = pkin(3) * t273;
t328 = r_i_i_C(3) * t269;
t325 = t269 * t275;
t324 = t269 * t277;
t323 = t275 * t268;
t321 = t275 * t283;
t320 = qJD(1) * t280;
t319 = qJD(1) * t283;
t311 = t268 * t321;
t310 = t268 * t319;
t309 = t268 * t320;
t298 = t343 * t268;
t296 = t302 * t278;
t294 = -t268 * t329 - t328;
t293 = -r_i_i_C(2) * t281 - t278 * t338;
t292 = -r_i_i_C(3) * t268 - t269 * t343;
t291 = -t298 - t324;
t290 = t277 * t312 + t344 * t280 + t310 * t329 + t319 * t328;
t289 = t277 * t311 + t344 * t283 + t343 * t309 + t320 * t324;
t288 = t301 * t280 + t311;
t282 = cos(qJ(2));
t287 = -t273 * t331 + t292 * t275 - t282 * t326;
t286 = t275 * (t292 - t332);
t285 = pkin(5) * t317 + qJD(4) + (-pkin(2) * t282 - t327 * t268 - t269 * t270 - pkin(1) - t332) * qJD(1);
t240 = t288 * t278 - t283 * t297;
t284 = r_i_i_C(3) * t325 + (-t298 - t333) * t275 + (t293 * qJD(5) - t275 * t277) * t269 + t345 * t268;
t264 = -pkin(2) * t279 - t333;
t243 = -t281 * t295 + (t281 * t323 + t296) * t280;
t241 = t288 * t281 + t283 * t296;
t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) - t339 * t280 + t285 * t283 (-t264 + t294) * t320 + t287 * t283 + t289 (t294 + t333) * t320 + t283 * t286 + t289, t319, t241 * r_i_i_C(2) + t338 * t240, t269 * t321 - t309; -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t285 * t280 + t339 * t283, t287 * t280 + (t264 + t291) * t319 + t290, t280 * t286 + (t291 - t333) * t319 + t290, t320, t243 * r_i_i_C(2) - t338 * t242, t269 * t322 + t310; 0, t284 - t313, t284, 0, t293 * t325 + (-t281 * t338 + t329) * t268 * qJD(5), t323;];
JaD_transl  = t1;
