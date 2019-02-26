% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:51
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR8_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR8_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR8_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RRRRRR8_jacobiaD_transl_3_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:51:24
% EndTime: 2019-02-26 22:51:24
% DurationCPUTime: 0.35s
% Computational Cost: add. (229->67), mult. (725->118), div. (0->0), fcn. (712->10), ass. (0->53)
t302 = cos(pkin(7));
t304 = sin(qJ(3));
t307 = cos(qJ(3));
t321 = r_i_i_C(1) * t304 + r_i_i_C(2) * t307;
t300 = sin(pkin(7));
t340 = pkin(10) + r_i_i_C(3);
t342 = t340 * t300;
t345 = -t302 * t321 + t342;
t303 = cos(pkin(6));
t305 = sin(qJ(2));
t309 = cos(qJ(1));
t329 = t309 * t305;
t306 = sin(qJ(1));
t308 = cos(qJ(2));
t331 = t306 * t308;
t292 = t303 * t329 + t331;
t317 = t303 * t331 + t329;
t289 = qJD(1) * t317 + qJD(2) * t292;
t301 = sin(pkin(6));
t323 = qJD(1) * t300 * t301;
t313 = -qJD(3) * t292 - t289 * t302 + t306 * t323;
t344 = r_i_i_C(1) * t313;
t341 = t302 * t340 + pkin(9);
t339 = t301 * t306;
t338 = t301 * t309;
t337 = t302 * t304;
t336 = t302 * t307;
t335 = t304 * t305;
t334 = t304 * t308;
t333 = t305 * t307;
t332 = t306 * t305;
t330 = t307 * t308;
t328 = t309 * t308;
t327 = qJD(3) * t300;
t324 = t303 * t332;
t322 = t327 * t338;
t290 = -qJD(1) * t324 - qJD(2) * t332 + (qJD(2) * t303 + qJD(1)) * t328;
t291 = -t303 * t328 + t332;
t320 = qJD(3) * t291 * t302 - t290;
t319 = t307 * r_i_i_C(1) - t304 * r_i_i_C(2) + pkin(2);
t318 = t300 * t339 - t302 * t317;
t316 = t324 - t328;
t287 = qJD(1) * t291 + qJD(2) * t316;
t315 = t287 * t302 + t309 * t323;
t314 = t320 + t322;
t312 = t313 * r_i_i_C(2);
t311 = (-t302 * t333 - t334) * r_i_i_C(1) + (t302 * t335 - t330) * r_i_i_C(2);
t310 = (-t302 * t334 - t333) * r_i_i_C(1) + (-t302 * t330 + t335) * r_i_i_C(2);
t295 = t307 * t322;
t288 = qJD(1) * t292 + qJD(2) * t317;
t286 = -t288 * t307 + t315 * t304 + (t304 * t316 + t307 * t318) * qJD(3);
t285 = t288 * t304 + t315 * t307 + (-t304 * t318 + t307 * t316) * qJD(3);
t1 = [-t290 * pkin(2) + t295 * r_i_i_C(1) - t289 * t342 + (-t309 * pkin(1) - t341 * t339) * qJD(1) + (r_i_i_C(1) * t320 - t312) * t307 + (-r_i_i_C(2) * t314 - t344) * t304, t319 * t287 - t345 * t288 + ((t304 * t317 + t316 * t336) * r_i_i_C(1) + (t307 * t317 - t316 * t337) * r_i_i_C(2)) * qJD(3), t285 * r_i_i_C(1) - t286 * r_i_i_C(2), 0, 0, 0; -t288 * pkin(2) + t286 * r_i_i_C(1) + t285 * r_i_i_C(2) - t287 * t342 + (-t306 * pkin(1) + t341 * t338) * qJD(1), -t319 * t289 + t345 * t290 + ((t291 * t304 - t292 * t336) * r_i_i_C(1) + (t291 * t307 + t292 * t337) * r_i_i_C(2)) * qJD(3), t295 * r_i_i_C(2) + (r_i_i_C(2) * t320 + t344) * t307 + (r_i_i_C(1) * t314 - t312) * t304, 0, 0, 0; 0 (t311 * qJD(3) + (-t305 * pkin(2) + t308 * t342 + t310) * qJD(2)) * t301, -t321 * t303 * t327 + (qJD(2) * t311 + qJD(3) * t310) * t301, 0, 0, 0;];
JaD_transl  = t1;
