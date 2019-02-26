% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR6_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:06:34
% EndTime: 2019-02-26 22:06:34
% DurationCPUTime: 0.36s
% Computational Cost: add. (459->86), mult. (980->134), div. (0->0), fcn. (944->10), ass. (0->53)
t301 = qJ(3) + pkin(11);
t299 = sin(t301);
t300 = cos(t301);
t305 = sin(qJ(3));
t339 = r_i_i_C(3) + qJ(5);
t342 = pkin(4) - r_i_i_C(2);
t347 = (pkin(3) * t305 + t342 * t299 - t339 * t300) * qJD(3) - t299 * qJD(5);
t307 = sin(qJ(1));
t303 = cos(pkin(6));
t320 = qJD(2) * t303 + qJD(1);
t306 = sin(qJ(2));
t333 = t307 * t306;
t325 = t303 * t333;
t328 = qJD(2) * t306;
t309 = cos(qJ(2));
t310 = cos(qJ(1));
t330 = t310 * t309;
t284 = -qJD(1) * t325 - t307 * t328 + t320 * t330;
t331 = t310 * t306;
t332 = t307 * t309;
t289 = t303 * t331 + t332;
t302 = sin(pkin(6));
t334 = t302 * t310;
t319 = t289 * t299 + t300 * t334;
t329 = qJD(1) * t302;
t323 = t307 * t329;
t346 = t319 * qJD(3) - t284 * t300 - t299 * t323;
t318 = -t289 * t300 + t299 * t334;
t345 = t318 * qJD(3) - t284 * t299 + t300 * t323;
t308 = cos(qJ(3));
t298 = t308 * pkin(3) + pkin(2);
t344 = t339 * t299 + t342 * t300 + t298;
t340 = r_i_i_C(1) + qJ(4) + pkin(9);
t291 = -t325 + t330;
t338 = t291 * t299;
t337 = t302 * t306;
t336 = t302 * t307;
t335 = t302 * t308;
t327 = qJD(2) * t309;
t324 = t303 * t330;
t322 = t310 * t329;
t321 = t302 * t327;
t317 = t291 * t300 + t299 * t336;
t316 = t303 * t299 + t300 * t337;
t315 = t303 * t332 + t331;
t288 = -t324 + t333;
t285 = t316 * qJD(3) + t299 * t321;
t283 = t315 * qJD(1) + t289 * qJD(2);
t282 = t289 * qJD(1) + t315 * qJD(2);
t281 = -qJD(1) * t324 - t310 * t327 + t320 * t333;
t276 = t299 * t322 - qJD(3) * t338 + (qJD(3) * t336 - t282) * t300;
t275 = t317 * qJD(3) - t282 * t299 - t300 * t322;
t1 = [-t319 * qJD(5) - t284 * t298 - t288 * qJD(4) - t340 * t283 + t342 * t346 + t339 * t345 + (-t310 * pkin(1) - pkin(8) * t336) * qJD(1) + (-t305 * t323 + (t289 * t305 + t308 * t334) * qJD(3)) * pkin(3), t291 * qJD(4) + t344 * t281 - t340 * t282 + t315 * t347, t317 * qJD(5) + t339 * t276 - t342 * t275 + (t308 * t322 + t282 * t305 + (-t291 * t308 - t305 * t336) * qJD(3)) * pkin(3), -t281, t275, 0; -(t300 * t336 - t338) * qJD(5) - t282 * t298 + t315 * qJD(4) - t340 * t281 + t342 * t276 + t339 * t275 + (-t307 * pkin(1) + pkin(8) * t334) * qJD(1) + (t305 * t322 + (-t291 * t305 + t307 * t335) * qJD(3)) * pkin(3), t289 * qJD(4) - t283 * t344 + t340 * t284 + t347 * t288, -t318 * qJD(5) - t339 * t346 + t342 * t345 + (t308 * t323 - t284 * t305 + (-t289 * t308 + t305 * t334) * qJD(3)) * pkin(3), t283, -t345, 0; 0 ((-qJD(2) * t344 + qJD(4)) * t306 + (t340 * qJD(2) - t347) * t309) * t302, t316 * qJD(5) + t339 * (t300 * t321 + (-t299 * t337 + t300 * t303) * qJD(3)) - t342 * t285 + (-t305 * t321 + (-t303 * t305 - t306 * t335) * qJD(3)) * pkin(3), t302 * t328, t285, 0;];
JaD_transl  = t1;
