% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:48
% EndTime: 2019-02-26 21:41:49
% DurationCPUTime: 0.53s
% Computational Cost: add. (753->95), mult. (1191->147), div. (0->0), fcn. (1113->10), ass. (0->63)
t302 = sin(qJ(2));
t295 = cos(pkin(10)) * pkin(3) + pkin(2);
t298 = pkin(10) + qJ(4);
t296 = sin(t298);
t297 = cos(t298);
t304 = cos(qJ(6));
t301 = sin(qJ(6));
t346 = -pkin(4) - pkin(5);
t323 = t301 * r_i_i_C(2) + t346;
t312 = t304 * r_i_i_C(1) - t323;
t324 = -t301 * r_i_i_C(1) - qJ(5);
t313 = t304 * r_i_i_C(2) - t324;
t347 = t313 * t296 + t312 * t297 + t295;
t305 = cos(qJ(2));
t331 = -r_i_i_C(3) - pkin(9) + pkin(8) + qJ(3);
t353 = t331 * t305;
t358 = t347 * t302 - t353;
t357 = (-t302 * t295 + t353) * qJD(2) + t302 * qJD(3);
t306 = cos(qJ(1));
t342 = t306 * t297;
t303 = sin(qJ(1));
t344 = t303 * t305;
t279 = t296 * t344 + t342;
t343 = t306 * t296;
t280 = t297 * t344 - t343;
t318 = t279 * t301 + t280 * t304;
t319 = t279 * t304 - t280 * t301;
t355 = (t318 * r_i_i_C(1) + t319 * r_i_i_C(2)) * qJD(6);
t314 = t296 * t301 + t297 * t304;
t315 = t296 * t304 - t297 * t301;
t354 = (t312 * t296 - t313 * t297) * qJD(4) - (t315 * r_i_i_C(1) - t314 * r_i_i_C(2)) * qJD(6) - t331 * qJD(2) - t296 * qJD(5);
t351 = (-qJD(4) + qJD(6)) * t302;
t335 = qJD(4) * t303;
t340 = qJD(1) * t306;
t311 = t296 * t340 + t297 * t335;
t334 = qJD(4) * t306;
t326 = t296 * t334;
t339 = qJD(2) * t302;
t329 = t303 * t339;
t341 = qJD(1) * t303;
t277 = -t296 * t329 - t297 * t341 + t311 * t305 - t326;
t274 = t277 * t304;
t338 = qJD(2) * t305;
t337 = qJD(2) * t306;
t336 = qJD(4) * t302;
t330 = pkin(3) * sin(pkin(10)) + pkin(7);
t328 = t296 * t335;
t327 = t302 * t337;
t325 = t297 * t334;
t320 = -(t314 * t351 - t315 * t338) * r_i_i_C(1) - (t314 * t338 + t315 * t351) * r_i_i_C(2);
t281 = -t303 * t297 + t305 * t343;
t282 = t303 * t296 + t305 * t342;
t317 = t281 * t304 - t282 * t301;
t316 = t281 * t301 + t282 * t304;
t310 = -t295 * t305 - t331 * t302 - pkin(1);
t308 = -qJD(2) * t347 + qJD(3);
t307 = t354 * t302 + t308 * t305;
t278 = t282 * qJD(1) - t297 * t329 - t305 * t328 - t325;
t276 = t305 * t326 + (t305 * t341 + t327) * t297 - t311;
t275 = t279 * qJD(1) + t296 * t327 - t305 * t325 - t328;
t271 = t317 * qJD(6) - t275 * t301 - t276 * t304;
t270 = -t316 * qJD(6) - t275 * t304 + t276 * t301;
t1 = [-t274 * r_i_i_C(2) - t279 * qJD(5) + t324 * t277 - t312 * t278 + (-t319 * r_i_i_C(1) + t318 * r_i_i_C(2)) * qJD(6) - t357 * t303 + (-t330 * t303 + t310 * t306) * qJD(1), t307 * t306 + t358 * t341, -t302 * t341 + t305 * t337, t282 * qJD(5) - t313 * t276 + t312 * t275 + (t316 * r_i_i_C(1) + t317 * r_i_i_C(2)) * qJD(6), -t275, t270 * r_i_i_C(1) - t271 * r_i_i_C(2); t271 * r_i_i_C(1) + t270 * r_i_i_C(2) - t275 * qJ(5) + t281 * qJD(5) + t346 * t276 + t357 * t306 + (t310 * t303 + t330 * t306) * qJD(1), t307 * t303 - t358 * t340, t302 * t340 + t303 * t338, -t274 * r_i_i_C(1) + t280 * qJD(5) + t323 * t277 + t313 * t278 + t355, t277 (-t278 * t301 + t274) * r_i_i_C(1) + (-t277 * t301 - t278 * t304) * r_i_i_C(2) - t355; 0, t308 * t302 - t354 * t305, t339 (-qJ(5) * t336 + t346 * t338) * t296 + (qJ(5) * t338 + (t346 * qJD(4) + qJD(5)) * t302) * t297 - t320, t296 * t338 + t297 * t336, t320;];
JaD_transl  = t1;
