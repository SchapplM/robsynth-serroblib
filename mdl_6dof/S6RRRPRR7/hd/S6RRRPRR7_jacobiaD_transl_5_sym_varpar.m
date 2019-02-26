% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRR7
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:19
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRPRR7_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:19:19
% EndTime: 2019-02-26 22:19:20
% DurationCPUTime: 0.31s
% Computational Cost: add. (474->78), mult. (741->117), div. (0->0), fcn. (687->12), ass. (0->57)
t299 = sin(qJ(1));
t296 = cos(pkin(6));
t311 = qJD(2) * t296 + qJD(1);
t298 = sin(qJ(2));
t327 = t299 * t298;
t314 = t296 * t327;
t318 = qJD(2) * t298;
t301 = cos(qJ(2));
t302 = cos(qJ(1));
t324 = t302 * t301;
t269 = -qJD(1) * t314 - t299 * t318 + t311 * t324;
t293 = qJD(3) + qJD(5);
t295 = sin(pkin(6));
t328 = t295 * t302;
t335 = t293 * t328 - t269;
t294 = qJ(3) + pkin(12);
t279 = pkin(4) * cos(t294) + cos(qJ(3)) * pkin(3);
t278 = sin(qJ(3)) * pkin(3) + pkin(4) * sin(t294);
t275 = t278 * qJD(3);
t290 = qJ(5) + t294;
t286 = sin(t290);
t287 = cos(t290);
t310 = r_i_i_C(1) * t286 + r_i_i_C(2) * t287;
t334 = t310 * t293 + t275;
t333 = pkin(8) + t278;
t332 = r_i_i_C(3) + pkin(10) + qJ(4) + pkin(9);
t325 = t302 * t298;
t326 = t299 * t301;
t272 = t296 * t325 + t326;
t331 = t272 * t293;
t330 = t286 * t293;
t329 = t295 * t299;
t274 = -t314 + t324;
t319 = qJD(1) * t302;
t305 = -t274 * t293 + t295 * t319;
t307 = t296 * t326 + t325;
t267 = t272 * qJD(1) + t307 * qJD(2);
t309 = t293 * t329 - t267;
t262 = -t309 * t286 + t305 * t287;
t263 = t305 * t286 + t309 * t287;
t323 = t262 * r_i_i_C(1) - t263 * r_i_i_C(2);
t320 = qJD(1) * t299;
t306 = t295 * t320 - t331;
t312 = t335 * t287;
t322 = (t335 * t286 + t306 * t287) * r_i_i_C(1) + (-t306 * t286 + t312) * r_i_i_C(2);
t317 = qJD(2) * t301;
t304 = -t293 * t296 - t295 * t317;
t316 = t293 * t295 * t298;
t321 = (t304 * t286 - t287 * t316) * r_i_i_C(1) + (t286 * t316 + t304 * t287) * r_i_i_C(2);
t313 = t296 * t324;
t277 = pkin(2) + t279;
t308 = t287 * r_i_i_C(1) - t286 * r_i_i_C(2) + t277;
t276 = t279 * qJD(3);
t271 = -t313 + t327;
t268 = t307 * qJD(1) + t272 * qJD(2);
t266 = -qJD(1) * t313 - t302 * t317 + t311 * t327;
t1 = [(t272 * t330 + t312) * r_i_i_C(1) + (t269 * t286 + t287 * t331) * r_i_i_C(2) - t269 * t277 + t272 * t275 - t271 * qJD(4) - pkin(1) * t319 - t332 * t268 + ((-r_i_i_C(2) * t330 + t276) * t302 + (-t310 - t333) * t320) * t295, t274 * qJD(4) + t308 * t266 - t332 * t267 + t307 * t334, t267 * t278 - t274 * t276 + (-t275 * t299 + t279 * t319) * t295 + t323, -t266, t323, 0; t276 * t329 + t263 * r_i_i_C(1) + t262 * r_i_i_C(2) + t307 * qJD(4) - t267 * t277 - t274 * t275 - t332 * t266 + (-pkin(1) * t299 + t333 * t328) * qJD(1), t272 * qJD(4) - t308 * t268 + t332 * t269 + t334 * t271, -t269 * t278 - t272 * t276 + (t275 * t302 + t279 * t320) * t295 + t322, t268, t322, 0; 0 (qJD(4) * t298 - t334 * t301 + (-t308 * t298 + t332 * t301) * qJD(2)) * t295, -t296 * t275 + (-t276 * t298 - t278 * t317) * t295 + t321, t295 * t318, t321, 0;];
JaD_transl  = t1;
