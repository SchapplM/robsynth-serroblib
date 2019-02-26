% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:03
% EndTime: 2019-02-26 21:02:03
% DurationCPUTime: 0.49s
% Computational Cost: add. (712->91), mult. (1162->142), div. (0->0), fcn. (1084->10), ass. (0->58)
t298 = sin(qJ(3));
t301 = cos(qJ(3));
t324 = -r_i_i_C(3) - pkin(9) + pkin(8);
t343 = t324 * t301;
t349 = (-t298 * pkin(3) + t343) * qJD(3);
t297 = sin(qJ(4));
t300 = cos(qJ(4));
t299 = cos(qJ(6));
t296 = sin(qJ(6));
t336 = -pkin(4) - pkin(5);
t316 = t296 * r_i_i_C(2) + t336;
t307 = t299 * r_i_i_C(1) - t316;
t319 = -t296 * r_i_i_C(1) - qJ(5);
t308 = t299 * r_i_i_C(2) - t319;
t337 = t308 * t297 + t307 * t300 + pkin(3);
t348 = t337 * t298 - t343;
t347 = (qJD(4) - qJD(6)) * t298;
t295 = qJ(1) + pkin(10);
t293 = sin(t295);
t294 = cos(t295);
t333 = t297 * t301;
t277 = t293 * t333 + t294 * t300;
t332 = t300 * t301;
t334 = t294 * t297;
t278 = t293 * t332 - t334;
t313 = t277 * t296 + t278 * t299;
t314 = t277 * t299 - t278 * t296;
t345 = (r_i_i_C(1) * t313 + r_i_i_C(2) * t314) * qJD(6);
t309 = t296 * t297 + t299 * t300;
t310 = t296 * t300 - t297 * t299;
t344 = -t297 * qJD(5) - t324 * qJD(3) + (t310 * r_i_i_C(1) + t309 * r_i_i_C(2)) * qJD(6) + (t307 * t297 - t308 * t300) * qJD(4);
t326 = qJD(4) * t300;
t330 = qJD(1) * t294;
t305 = t293 * t326 + t297 * t330;
t329 = qJD(3) * t298;
t322 = t293 * t329;
t331 = qJD(1) * t293;
t323 = t300 * t331;
t327 = qJD(4) * t297;
t275 = -t294 * t327 - t297 * t322 + t305 * t301 - t323;
t272 = t275 * t299;
t328 = qJD(3) * t301;
t321 = t294 * t326;
t320 = t301 * t327;
t315 = -(-t309 * t347 + t310 * t328) * r_i_i_C(1) - (t309 * t328 + t310 * t347) * r_i_i_C(2);
t279 = -t293 * t300 + t294 * t333;
t280 = t293 * t297 + t294 * t332;
t312 = t279 * t299 - t280 * t296;
t311 = t279 * t296 + t280 * t299;
t306 = -pkin(3) * t301 - t324 * t298 - pkin(2);
t303 = qJD(3) * t337;
t302 = t344 * t298 - t301 * t303;
t276 = t280 * qJD(1) - t293 * t320 - t300 * t322 - t321;
t274 = t301 * t323 + (t300 * t329 + t320) * t294 - t305;
t273 = t277 * qJD(1) - t293 * t327 - t301 * t321 + t329 * t334;
t269 = t312 * qJD(6) - t273 * t296 - t274 * t299;
t268 = -t311 * qJD(6) - t273 * t299 + t274 * t296;
t1 = [-t272 * r_i_i_C(2) - t277 * qJD(5) + t319 * t275 - t307 * t276 + (-t314 * r_i_i_C(1) + t313 * r_i_i_C(2)) * qJD(6) - t293 * t349 + (-cos(qJ(1)) * pkin(1) - t293 * pkin(7) + t306 * t294) * qJD(1), 0, t302 * t294 + t348 * t331, t280 * qJD(5) - t308 * t274 + t307 * t273 + (r_i_i_C(1) * t311 + r_i_i_C(2) * t312) * qJD(6), -t273, t268 * r_i_i_C(1) - t269 * r_i_i_C(2); t269 * r_i_i_C(1) + t268 * r_i_i_C(2) - t273 * qJ(5) + t279 * qJD(5) + t336 * t274 + t294 * t349 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t294 + t306 * t293) * qJD(1), 0, t302 * t293 - t348 * t330, -t272 * r_i_i_C(1) + t278 * qJD(5) + t316 * t275 + t308 * t276 + t345, t275 (-t276 * t296 + t272) * r_i_i_C(1) + (-t275 * t296 - t276 * t299) * r_i_i_C(2) - t345; 0, 0, -t298 * t303 - t344 * t301 (t300 * qJ(5) + t336 * t297) * t328 + (qJD(5) * t300 + (-t297 * qJ(5) + t336 * t300) * qJD(4)) * t298 - t315, t297 * t328 + t298 * t326, t315;];
JaD_transl  = t1;
