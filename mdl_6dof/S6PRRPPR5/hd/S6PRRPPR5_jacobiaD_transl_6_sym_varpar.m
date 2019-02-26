% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d6,theta1,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:00:38
% EndTime: 2019-02-26 20:00:39
% DurationCPUTime: 0.39s
% Computational Cost: add. (399->71), mult. (1032->125), div. (0->0), fcn. (1036->12), ass. (0->51)
t311 = sin(qJ(3));
t313 = cos(qJ(3));
t306 = pkin(11) + qJ(6);
t304 = sin(t306);
t305 = cos(t306);
t328 = t305 * r_i_i_C(1) - t304 * r_i_i_C(2);
t317 = qJD(6) * t328 + qJD(4);
t327 = -t304 * r_i_i_C(1) - t305 * r_i_i_C(2);
t318 = pkin(5) * sin(pkin(11)) + qJ(4) - t327;
t336 = pkin(3) + r_i_i_C(3) + pkin(9) + qJ(5);
t315 = (t311 * t336 - t313 * t318) * qJD(3) - t313 * qJD(5) - t317 * t311;
t308 = sin(pkin(10));
t312 = sin(qJ(2));
t314 = cos(qJ(2));
t343 = cos(pkin(10));
t344 = cos(pkin(6));
t326 = t344 * t343;
t294 = t308 * t314 + t312 * t326;
t309 = sin(pkin(6));
t332 = t309 * t343;
t346 = t294 * t313 - t311 * t332;
t333 = t308 * t344;
t296 = -t312 * t333 + t343 * t314;
t341 = t309 * t311;
t340 = t309 * t313;
t339 = t309 * t314;
t338 = qJD(2) * t312;
t335 = t309 * t338;
t334 = qJD(2) * t339;
t325 = t314 * t326;
t324 = -t296 * t311 + t308 * t340;
t323 = t296 * t313 + t308 * t341;
t322 = t327 * qJD(6);
t321 = pkin(8) + cos(pkin(11)) * pkin(5) + pkin(4) + t328;
t297 = t312 * t341 - t313 * t344;
t320 = t311 * t344 + t312 * t340;
t319 = -t294 * t311 - t313 * t332;
t295 = t312 * t343 + t314 * t333;
t316 = -t311 * t318 - t313 * t336 - pkin(2);
t293 = t308 * t312 - t325;
t292 = t296 * qJD(2);
t291 = t295 * qJD(2);
t290 = t294 * qJD(2);
t289 = -qJD(2) * t325 + t308 * t338;
t288 = -qJD(3) * t297 + t313 * t334;
t287 = qJD(3) * t320 + t311 * t334;
t282 = qJD(3) * t324 - t291 * t313;
t281 = qJD(3) * t323 - t291 * t311;
t280 = qJD(3) * t319 - t289 * t313;
t279 = t346 * qJD(3) - t289 * t311;
t1 = [0, -t291 * t321 + t292 * t316 + t295 * t315 + t296 * t322, qJD(5) * t324 - t281 * t336 + t282 * t318 + t317 * t323, t281, t282 (t281 * t305 - t292 * t304) * r_i_i_C(1) + (-t281 * t304 - t292 * t305) * r_i_i_C(2) + ((-t295 * t305 + t304 * t324) * r_i_i_C(1) + (t295 * t304 + t305 * t324) * r_i_i_C(2)) * qJD(6); 0, -t289 * t321 + t290 * t316 + t293 * t315 + t294 * t322, qJD(5) * t319 - t336 * t279 + t318 * t280 + t317 * t346, t279, t280 (t279 * t305 - t290 * t304) * r_i_i_C(1) + (-t279 * t304 - t290 * t305) * r_i_i_C(2) + ((-t293 * t305 + t304 * t319) * r_i_i_C(1) + (t293 * t304 + t305 * t319) * r_i_i_C(2)) * qJD(6); 0 ((qJD(2) * t316 + t322) * t312 + (t321 * qJD(2) - t315) * t314) * t309, -t297 * qJD(5) - t287 * t336 + t288 * t318 + t317 * t320, t287, t288 (t287 * t305 - t304 * t335) * r_i_i_C(1) + (-t287 * t304 - t305 * t335) * r_i_i_C(2) + ((-t297 * t304 + t305 * t339) * r_i_i_C(1) + (-t297 * t305 - t304 * t339) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
