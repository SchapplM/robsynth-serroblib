% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRRPP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRRPP1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRPP1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRRPP1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPP1_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:08:46
% EndTime: 2019-02-26 20:08:46
% DurationCPUTime: 0.48s
% Computational Cost: add. (385->92), mult. (978->164), div. (0->0), fcn. (980->12), ass. (0->56)
t306 = sin(qJ(3));
t309 = cos(qJ(3));
t301 = qJ(4) + pkin(11);
t299 = sin(t301);
t300 = cos(t301);
t308 = cos(qJ(4));
t348 = t308 * pkin(4) + t300 * r_i_i_C(1) - t299 * r_i_i_C(2);
t319 = pkin(3) + t348;
t342 = r_i_i_C(3) + qJ(5) + pkin(9);
t349 = -(t319 * t306 - t342 * t309) * qJD(3) + t306 * qJD(5);
t302 = sin(pkin(10));
t307 = sin(qJ(2));
t310 = cos(qJ(2));
t340 = cos(pkin(10));
t341 = cos(pkin(6));
t320 = t341 * t340;
t289 = t302 * t310 + t307 * t320;
t303 = sin(pkin(6));
t326 = t303 * t340;
t279 = t289 * t309 - t306 * t326;
t327 = t302 * t341;
t291 = -t307 * t327 + t340 * t310;
t305 = sin(qJ(4));
t316 = t305 * pkin(4) + t299 * r_i_i_C(1) + t300 * r_i_i_C(2);
t338 = t303 * t306;
t337 = t303 * t309;
t336 = t303 * t310;
t335 = qJD(2) * t307;
t334 = qJD(4) * t299;
t333 = qJD(4) * t300;
t332 = qJD(4) * t308;
t331 = qJD(4) * t309;
t329 = t303 * t335;
t328 = qJD(2) * t336;
t318 = t310 * t320;
t317 = -t291 * t306 + t302 * t337;
t281 = t291 * t309 + t302 * t338;
t315 = -t289 * t306 - t309 * t326;
t314 = -t307 * t338 + t341 * t309;
t293 = t341 * t306 + t307 * t337;
t313 = qJD(4) * t316;
t290 = t340 * t307 + t310 * t327;
t312 = -t342 * t306 - t319 * t309 - pkin(2);
t311 = t316 * t331 - t349;
t288 = t302 * t307 - t318;
t287 = t291 * qJD(2);
t286 = t290 * qJD(2);
t285 = t289 * qJD(2);
t284 = -qJD(2) * t318 + t302 * t335;
t283 = t314 * qJD(3) + t309 * t328;
t282 = t293 * qJD(3) + t306 * t328;
t277 = t317 * qJD(3) - t286 * t309;
t276 = t281 * qJD(3) - t286 * t306;
t275 = t315 * qJD(3) - t284 * t309;
t274 = t279 * qJD(3) - t284 * t306;
t1 = [0 (-t286 * t299 + t291 * t333) * r_i_i_C(1) + (-t286 * t300 - t291 * t334) * r_i_i_C(2) - t286 * pkin(8) + (-t286 * t305 + t291 * t332) * pkin(4) + t312 * t287 + t311 * t290, t281 * qJD(5) - t319 * t276 + t342 * t277 - t317 * t313 (-t277 * t299 + t287 * t300) * r_i_i_C(1) + (-t277 * t300 - t287 * t299) * r_i_i_C(2) + ((-t281 * t300 - t290 * t299) * r_i_i_C(1) + (t281 * t299 - t290 * t300) * r_i_i_C(2)) * qJD(4) + (-t277 * t305 + t287 * t308 + (-t281 * t308 - t290 * t305) * qJD(4)) * pkin(4), t276, 0; 0 (-t284 * t299 + t289 * t333) * r_i_i_C(1) + (-t284 * t300 - t289 * t334) * r_i_i_C(2) - t284 * pkin(8) + (-t284 * t305 + t289 * t332) * pkin(4) + t312 * t285 + t311 * t288, t279 * qJD(5) - t319 * t274 + t342 * t275 - t315 * t313 (-t275 * t299 + t285 * t300) * r_i_i_C(1) + (-t275 * t300 - t285 * t299) * r_i_i_C(2) + ((-t279 * t300 - t288 * t299) * r_i_i_C(1) + (t279 * t299 - t288 * t300) * r_i_i_C(2)) * qJD(4) + (-t275 * t305 + t285 * t308 + (-t279 * t308 - t288 * t305) * qJD(4)) * pkin(4), t274, 0; 0 ((t312 * qJD(2) + t348 * qJD(4)) * t307 + (qJD(2) * pkin(8) + t316 * (qJD(2) - t331) + t349) * t310) * t303, t293 * qJD(5) - t319 * t282 + t342 * t283 - t314 * t313 (-t283 * t299 + t300 * t329) * r_i_i_C(1) + (-t283 * t300 - t299 * t329) * r_i_i_C(2) + ((-t293 * t300 + t299 * t336) * r_i_i_C(1) + (t293 * t299 + t300 * t336) * r_i_i_C(2)) * qJD(4) + (t308 * t329 - t283 * t305 + (-t293 * t308 + t305 * t336) * qJD(4)) * pkin(4), t282, 0;];
JaD_transl  = t1;
