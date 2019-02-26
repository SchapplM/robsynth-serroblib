% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRPRR9
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR9_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR9_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR9_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRPRR9_jacobiaD_transl_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:53:25
% EndTime: 2019-02-26 20:53:25
% DurationCPUTime: 0.35s
% Computational Cost: add. (208->87), mult. (693->156), div. (0->0), fcn. (698->12), ass. (0->61)
t338 = pkin(9) + qJ(4);
t303 = sin(pkin(12));
t311 = sin(qJ(1));
t337 = t303 * t311;
t304 = sin(pkin(7));
t305 = sin(pkin(6));
t336 = t304 * t305;
t308 = cos(pkin(7));
t310 = sin(qJ(3));
t335 = t308 * t310;
t309 = cos(pkin(6));
t313 = cos(qJ(1));
t334 = t309 * t313;
t307 = cos(pkin(12));
t333 = t311 * t307;
t332 = qJD(1) * t311;
t331 = qJD(1) * t313;
t330 = qJD(3) * t310;
t312 = cos(qJ(3));
t329 = qJD(3) * t312;
t328 = pkin(3) * t330;
t327 = t309 * t337;
t326 = t305 * t332;
t325 = t305 * t331;
t324 = t304 * t330;
t323 = t304 * t329;
t322 = t308 * t329;
t302 = sin(pkin(13));
t306 = cos(pkin(13));
t319 = t302 * t312 + t306 * t310;
t275 = t319 * t304;
t321 = pkin(3) * t304 * t310 + r_i_i_C(1) * t275 + t338 * t308 + qJ(2);
t270 = t302 * t324 - t306 * t323;
t320 = pkin(3) * t323 - r_i_i_C(1) * t270 + qJD(4) * t308 + qJD(2);
t292 = t302 * t310 - t312 * t306;
t284 = -t307 * t334 + t337;
t318 = t303 * t313 + t309 * t333;
t285 = t303 * t334 + t333;
t317 = qJD(3) * t319;
t272 = t302 * t308 * t330 - t306 * t322;
t281 = -qJD(1) * t327 + t307 * t331;
t289 = -t302 * t329 - t306 * t330;
t316 = -t284 * t272 + t281 * t292 - t285 * t289;
t273 = t308 * t317;
t288 = t292 * qJD(3);
t315 = -t284 * t273 + t281 * t319 - t285 * t288;
t277 = t319 * t308;
t278 = t284 * qJD(1);
t279 = t285 * qJD(1);
t287 = t307 * t313 - t327;
t314 = -t272 * t318 - t278 * t277 - t279 * t292 - t287 * t289;
t301 = pkin(3) * t312 + pkin(2);
t291 = pkin(3) * t322 - qJD(4) * t304;
t283 = pkin(3) * t335 - t338 * t304;
t280 = t318 * qJD(1);
t276 = t292 * t308;
t274 = t292 * t304;
t271 = t304 * t317;
t269 = -t278 * t304 + t308 * t325;
t268 = t318 * t273 - t278 * t276 + t279 * t319 + t287 * t288 + (-t271 * t311 - t274 * t331) * t305;
t1 = [t316 * r_i_i_C(1) + t315 * r_i_i_C(2) - t281 * t301 + t285 * t328 + t284 * t291 - pkin(1) * t331 + (t277 * r_i_i_C(1) - t276 * r_i_i_C(2) - t304 * r_i_i_C(3) + t283) * t280 + ((-t271 * r_i_i_C(2) + t320) * t313 + (r_i_i_C(2) * t274 - r_i_i_C(3) * t308 - t321) * t332) * t305, t325, t268 * r_i_i_C(1) + ((t270 * t311 - t275 * t331) * t305 + t314) * r_i_i_C(2) + (t279 * t310 + (t278 * t308 + t304 * t325) * t312 + (-t287 * t312 + (t308 * t318 - t311 * t336) * t310) * qJD(3)) * pkin(3), t269, 0, 0; -t314 * r_i_i_C(1) + t268 * r_i_i_C(2) + t269 * r_i_i_C(3) - t279 * t301 - t287 * t328 + t278 * t283 - t318 * t291 - pkin(1) * t332 + (t320 * t311 + t321 * t331) * t305, t326 (t280 * t276 - t315) * r_i_i_C(1) + (t280 * t277 + t316) * r_i_i_C(2) + ((t271 * t313 - t274 * t332) * r_i_i_C(1) + (-t270 * t313 - t275 * t332) * r_i_i_C(2)) * t305 + (-t281 * t310 + (-t280 * t308 + t304 * t326) * t312 + (-t285 * t312 + (t284 * t308 + t313 * t336) * t310) * qJD(3)) * pkin(3), t280 * t304 + t308 * t326, 0, 0; 0, 0 (-pkin(3) * t324 - r_i_i_C(1) * t271 + r_i_i_C(2) * t270) * t309 + ((-t273 * t307 + t288 * t303) * r_i_i_C(1) + (t272 * t307 - t289 * t303) * r_i_i_C(2) + (-t303 * t312 - t307 * t335) * qJD(3) * pkin(3)) * t305, 0, 0, 0;];
JaD_transl  = t1;
