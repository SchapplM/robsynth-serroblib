% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRPP4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:27:01
% EndTime: 2019-02-26 22:27:01
% DurationCPUTime: 0.42s
% Computational Cost: add. (837->83), mult. (840->113), div. (0->0), fcn. (698->10), ass. (0->65)
t314 = qJ(3) + qJ(4);
t308 = pkin(10) + t314;
t306 = cos(t308);
t359 = r_i_i_C(3) + qJ(6);
t374 = t359 * t306;
t313 = qJD(3) + qJD(4);
t309 = sin(t314);
t315 = sin(qJ(3));
t358 = pkin(3) * qJD(3);
t362 = pkin(4) * t313;
t287 = -t309 * t362 - t315 * t358;
t305 = sin(t308);
t360 = r_i_i_C(2) + qJ(5) + pkin(9) + pkin(8);
t328 = t360 * qJD(2) + qJD(6) * t305 + t287;
t364 = pkin(5) + r_i_i_C(1);
t345 = t364 * t305;
t373 = (t345 - t374) * t313 - t328;
t310 = cos(t314);
t307 = pkin(4) * t310;
t318 = cos(qJ(3));
t301 = t318 * pkin(3) + t307;
t298 = pkin(2) + t301;
t363 = pkin(4) * t309;
t300 = t315 * pkin(3) + t363;
t316 = sin(qJ(2));
t319 = cos(qJ(2));
t372 = t328 * t319 + (pkin(7) + t300) * qJD(1) - (qJD(2) * t298 - qJD(5)) * t316;
t330 = -t359 * t305 - t364 * t306;
t325 = -t298 + t330;
t370 = t325 * t316 + t360 * t319;
t317 = sin(qJ(1));
t320 = cos(qJ(1));
t353 = t320 * t306;
t332 = t317 * t305 + t319 * t353;
t355 = t317 * t319;
t366 = t305 * t355 + t353;
t357 = t313 * t316;
t354 = t320 * t305;
t352 = qJD(1) * t317;
t351 = qJD(1) * t319;
t350 = qJD(1) * t320;
t349 = qJD(2) * t316;
t348 = qJD(2) * t319;
t347 = qJD(2) * t320;
t346 = t306 * qJD(6);
t343 = t313 * t354;
t341 = t316 * t346 + t348 * t374;
t339 = t317 * t349;
t338 = t316 * t347;
t337 = -t313 + t351;
t334 = t300 * t351 + t287;
t333 = t310 * (-t313 * t319 + qJD(1));
t331 = t317 * t313 * t306 + t305 * t350;
t281 = t366 * qJD(1) + t305 * t338 - t332 * t313;
t282 = t319 * t343 + (t317 * t351 + t338) * t306 - t331;
t327 = t332 * qJD(6) + t364 * t281 - t359 * t282;
t283 = -t305 * t339 - t306 * t352 + t331 * t319 - t343;
t284 = t332 * qJD(1) - t306 * t339 - t366 * t313;
t326 = -(-t306 * t355 + t354) * qJD(6) + t359 * t284 - t364 * t283;
t288 = t310 * t362 + t318 * t358;
t324 = qJD(1) * t301 - t288 * t319 + t300 * t349;
t323 = t325 * qJD(2) + qJD(5);
t322 = -t346 + t288 + (-t298 * t319 - t360 * t316 - pkin(1)) * qJD(1);
t321 = t316 * t373 + t323 * t319;
t1 = [-t359 * t283 - t364 * t284 - t317 * t372 + t322 * t320, t321 * t320 - t370 * t352, t334 * t317 + t324 * t320 + t327 (t320 * t333 + (t337 * t317 + t338) * t309) * pkin(4) + t327, -t316 * t352 + t319 * t347, -t281; -t359 * t281 - t364 * t282 + t322 * t317 + t320 * t372, t321 * t317 + t370 * t350, t324 * t317 - t334 * t320 + t326 (t317 * t333 + (-t337 * t320 + t339) * t309) * pkin(4) + t326, t316 * t350 + t317 * t348, t283; 0, t323 * t316 - t319 * t373 (-t300 - t345) * t348 + (t330 * t313 - t288) * t316 + t341 (-t345 - t363) * t348 + (t330 - t307) * t357 + t341, t349, t305 * t348 + t306 * t357;];
JaD_transl  = t1;
