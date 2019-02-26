% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:11
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:11:21
% EndTime: 2019-02-26 21:11:21
% DurationCPUTime: 0.34s
% Computational Cost: add. (708->72), mult. (744->103), div. (0->0), fcn. (626->9), ass. (0->57)
t313 = qJ(4) + qJ(5);
t310 = cos(t313);
t358 = r_i_i_C(3) + qJ(6);
t365 = t358 * t310;
t317 = cos(qJ(4));
t306 = t317 * pkin(4) + pkin(3);
t311 = pkin(10) + qJ(3);
t308 = cos(t311);
t309 = sin(t313);
t315 = sin(qJ(4));
t357 = pkin(4) * qJD(4);
t359 = r_i_i_C(2) + pkin(9) + pkin(8);
t322 = t359 * qJD(3) + qJD(6) * t309 - t315 * t357;
t307 = sin(t311);
t350 = qJD(3) * t307;
t364 = -t306 * t350 + t322 * t308;
t312 = qJD(4) + qJD(5);
t361 = pkin(5) + r_i_i_C(1);
t346 = t361 * t309;
t320 = (t346 - t365) * t312 - t322;
t329 = -t358 * t309 - t361 * t310;
t325 = -t306 + t329;
t360 = pkin(4) * t315;
t316 = sin(qJ(1));
t356 = t316 * t309;
t355 = t316 * t310;
t318 = cos(qJ(1));
t354 = t318 * t309;
t353 = t318 * t310;
t352 = qJD(1) * t316;
t351 = qJD(1) * t318;
t349 = qJD(3) * t308;
t348 = t310 * qJD(6);
t347 = t317 * t357;
t345 = t312 * t356;
t344 = t312 * t353;
t343 = t307 * t348 + t349 * t365;
t341 = t309 * t349;
t339 = t316 * t350;
t338 = t318 * t350;
t337 = pkin(7) + qJ(2) + t360;
t336 = qJD(1) * t308 - qJD(4);
t335 = (-qJD(4) * t308 + qJD(1)) * t317;
t333 = t308 * t353 + t356;
t332 = t310 * t352 + t312 * t354;
t331 = t309 * t351 + t312 * t355;
t330 = qJD(2) + t347 - t348;
t328 = -t306 * t308 - t359 * t307 - cos(pkin(10)) * pkin(2) - pkin(1);
t286 = t309 * t338 - t308 * t344 - t345 + (t308 * t356 + t353) * qJD(1);
t287 = t332 * t308 + t310 * t338 - t331;
t327 = t333 * qJD(6) + t361 * t286 - t358 * t287;
t288 = t331 * t308 - t309 * t339 - t332;
t289 = t333 * qJD(1) - t308 * t345 - t310 * t339 - t344;
t326 = -(-t308 * t355 + t354) * qJD(6) + t358 * t289 - t361 * t288;
t324 = t329 * t312;
t323 = qJD(3) * t325;
t1 = [-t361 * t289 - t358 * t288 + t330 * t318 - t364 * t316 + (-t337 * t316 + t328 * t318) * qJD(1), t351 (t318 * t323 - t359 * t352) * t308 + (t320 * t318 - t325 * t352) * t307 (t318 * t335 + (t336 * t316 + t338) * t315) * pkin(4) + t327, t327, -t286; -t361 * t287 - t358 * t286 + t330 * t316 + t364 * t318 + (t328 * t316 + t337 * t318) * qJD(1), t352 (t316 * t323 + t359 * t351) * t308 + (t320 * t316 + t325 * t351) * t307 (t316 * t335 + (-t336 * t318 + t339) * t315) * pkin(4) + t326, t326, t288; 0, 0, t307 * t323 - t320 * t308 (-t346 - t360) * t349 + (t324 - t347) * t307 + t343, t307 * t324 - t361 * t341 + t343, t307 * t312 * t310 + t341;];
JaD_transl  = t1;
