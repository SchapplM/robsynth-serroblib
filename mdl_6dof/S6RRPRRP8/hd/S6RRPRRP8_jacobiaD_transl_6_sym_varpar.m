% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:50
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:49:54
% EndTime: 2019-02-26 21:49:55
% DurationCPUTime: 0.38s
% Computational Cost: add. (777->73), mult. (767->104), div. (0->0), fcn. (649->10), ass. (0->57)
t304 = pkin(10) + qJ(4);
t302 = qJ(5) + t304;
t299 = cos(t302);
t348 = r_i_i_C(3) + qJ(6);
t362 = t348 * t299;
t305 = qJD(4) + qJD(5);
t298 = sin(t302);
t300 = sin(t304);
t347 = pkin(4) * qJD(4);
t349 = r_i_i_C(2) + pkin(9) + pkin(8) + qJ(3);
t314 = qJD(2) * t349 + qJD(6) * t298 - t300 * t347;
t352 = pkin(5) + r_i_i_C(1);
t333 = t352 * t298;
t361 = (t333 - t362) * t305 - t314;
t301 = cos(t304);
t291 = pkin(4) * t301 + cos(pkin(10)) * pkin(3) + pkin(2);
t306 = sin(qJ(2));
t308 = cos(qJ(2));
t351 = pkin(4) * t300;
t360 = t308 * t314 + qJD(1) * (pkin(7) + t351 + sin(pkin(10)) * pkin(3)) - (qJD(2) * t291 - qJD(3)) * t306;
t318 = -t348 * t298 - t299 * t352;
t315 = -t291 + t318;
t358 = t306 * t315 + t349 * t308;
t307 = sin(qJ(1));
t309 = cos(qJ(1));
t342 = t309 * t299;
t320 = t307 * t298 + t308 * t342;
t344 = t307 * t308;
t354 = t298 * t344 + t342;
t346 = t305 * t306;
t343 = t309 * t298;
t341 = qJD(1) * t307;
t340 = qJD(1) * t308;
t339 = qJD(1) * t309;
t338 = qJD(2) * t306;
t337 = qJD(2) * t308;
t336 = qJD(2) * t309;
t335 = t299 * qJD(6);
t334 = t301 * t347;
t332 = t299 * t346;
t330 = t305 * t343;
t328 = t306 * t335 + t337 * t362;
t326 = t307 * t338;
t325 = t306 * t336;
t324 = -qJD(4) + t340;
t322 = t301 * (-qJD(4) * t308 + qJD(1));
t319 = t307 * t305 * t299 + t298 * t339;
t277 = t354 * qJD(1) + t298 * t325 - t320 * t305;
t278 = t308 * t330 + (t307 * t340 + t325) * t299 - t319;
t317 = t320 * qJD(6) + t352 * t277 - t348 * t278;
t279 = -t298 * t326 - t299 * t341 + t308 * t319 - t330;
t280 = qJD(1) * t320 - t299 * t326 - t354 * t305;
t316 = -(-t299 * t344 + t343) * qJD(6) + t348 * t280 - t352 * t279;
t312 = qJD(2) * t315 + qJD(3);
t311 = t334 - t335 + (-t291 * t308 - t306 * t349 - pkin(1)) * qJD(1);
t310 = t361 * t306 + t312 * t308;
t1 = [-t348 * t279 - t352 * t280 - t360 * t307 + t311 * t309, t310 * t309 - t358 * t341, -t306 * t341 + t308 * t336 (t309 * t322 + (t307 * t324 + t325) * t300) * pkin(4) + t317, t317, -t277; -t348 * t277 - t352 * t278 + t311 * t307 + t360 * t309, t310 * t307 + t358 * t339, t306 * t339 + t307 * t337 (t307 * t322 + (-t324 * t309 + t326) * t300) * pkin(4) + t316, t316, t279; 0, t312 * t306 - t361 * t308, t338 (-t333 - t351) * t337 + (t305 * t318 - t334) * t306 + t328, -t352 * t332 + (-t337 * t352 - t348 * t346) * t298 + t328, t337 * t298 + t332;];
JaD_transl  = t1;
