% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPP6
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP6_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:28:16
% EndTime: 2019-02-26 22:28:16
% DurationCPUTime: 0.40s
% Computational Cost: add. (693->82), mult. (956->118), div. (0->0), fcn. (817->8), ass. (0->56)
t297 = qJ(3) + qJ(4);
t294 = sin(t297);
t295 = cos(t297);
t296 = qJD(3) + qJD(4);
t331 = pkin(4) + r_i_i_C(3) + qJ(6);
t317 = t331 * t294;
t298 = sin(qJ(3));
t344 = pkin(3) * qJD(3);
t329 = t298 * t344;
t330 = pkin(5) + r_i_i_C(1) + pkin(9) + pkin(8);
t345 = r_i_i_C(2) + qJ(5);
t305 = -t330 * qJD(2) - qJD(5) * t294 - qJD(6) * t295 + (-t345 * t295 + t317) * t296 + t329;
t301 = cos(qJ(3));
t293 = t301 * pkin(3) + pkin(2);
t316 = t331 * t295;
t311 = -t345 * t294 - t316;
t309 = -t293 + t311;
t299 = sin(qJ(2));
t302 = cos(qJ(2));
t351 = (-t293 * t299 + t330 * t302) * qJD(2) - t302 * t329;
t300 = sin(qJ(1));
t303 = cos(qJ(1));
t337 = t303 * t295;
t314 = t300 * t294 + t302 * t337;
t339 = t300 * t302;
t275 = t294 * t339 + t337;
t346 = pkin(3) * t298;
t342 = t296 * t299;
t340 = t300 * t295;
t338 = t303 * t294;
t336 = qJD(1) * t300;
t335 = qJD(1) * t302;
t334 = qJD(1) * t303;
t333 = qJD(2) * t299;
t332 = qJD(2) * t302;
t328 = t301 * t344;
t326 = t302 * t338;
t320 = t295 * t332;
t324 = t299 * t295 * qJD(5) + t345 * t320;
t323 = pkin(7) + t346;
t322 = t300 * t333;
t321 = t303 * t333;
t318 = -qJD(3) + t335;
t315 = (-qJD(3) * t302 + qJD(1)) * t301;
t313 = t294 * t334 + t296 * t340;
t312 = -t293 * t302 - t330 * t299 - pkin(1);
t308 = qJD(2) * t309;
t271 = t275 * qJD(1) + t294 * t321 - t314 * t296;
t272 = t296 * t326 + (t300 * t335 + t321) * t295 - t313;
t277 = t326 - t340;
t307 = t314 * qJD(5) - t277 * qJD(6) + t331 * t271 - t345 * t272;
t273 = -t294 * t322 - t295 * t336 - t296 * t338 + t313 * t302;
t274 = t314 * qJD(1) - t275 * t296 - t295 * t322;
t276 = t295 * t339 - t338;
t306 = t276 * qJD(5) - t275 * qJD(6) - t331 * t273 + t345 * t274;
t1 = [t303 * t328 - t275 * qJD(5) - t276 * qJD(6) - t345 * t273 - t331 * t274 - t351 * t300 + (-t323 * t300 + t312 * t303) * qJD(1) (t303 * t308 - t330 * t336) * t302 + (t305 * t303 - t309 * t336) * t299 (t303 * t315 + (t318 * t300 + t321) * t298) * pkin(3) + t307, t307, -t271, -t272; t300 * t328 + t277 * qJD(5) + t314 * qJD(6) - t345 * t271 - t331 * t272 + t351 * t303 + (t312 * t300 + t323 * t303) * qJD(1) (t300 * t308 + t330 * t334) * t302 + (t305 * t300 + t309 * t334) * t299 (t300 * t315 + (-t318 * t303 + t322) * t298) * pkin(3) + t306, t306, t273, t274; 0, t299 * t308 - t305 * t302 (-t317 - t346) * t332 + (-qJD(6) * t294 + t311 * t296 - t328) * t299 + t324, -t316 * t342 + ((-t345 * t296 - qJD(6)) * t299 - t331 * t332) * t294 + t324, t294 * t332 + t295 * t342, -t294 * t342 + t320;];
JaD_transl  = t1;
