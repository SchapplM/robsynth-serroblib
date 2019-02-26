% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:49
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:49:48
% EndTime: 2019-02-26 19:49:49
% DurationCPUTime: 0.33s
% Computational Cost: add. (287->69), mult. (917->127), div. (0->0), fcn. (916->10), ass. (0->51)
t297 = sin(qJ(6));
t300 = cos(qJ(6));
t312 = t297 * r_i_i_C(1) + t300 * r_i_i_C(2);
t308 = qJD(6) * t312;
t313 = t300 * r_i_i_C(1) - t297 * r_i_i_C(2);
t305 = qJD(6) * t313 + qJD(5);
t298 = sin(qJ(4));
t301 = cos(qJ(4));
t311 = qJ(5) + t312;
t322 = pkin(4) + pkin(9) + r_i_i_C(3);
t332 = t298 * t322 - t301 * t311 + qJ(3);
t331 = cos(pkin(6));
t294 = sin(pkin(10));
t296 = cos(pkin(10));
t299 = sin(qJ(2));
t302 = cos(qJ(2));
t317 = t302 * t331;
t284 = t294 * t317 + t296 * t299;
t330 = t284 * t301;
t295 = sin(pkin(6));
t329 = t295 * t298;
t328 = t295 * t299;
t327 = t295 * t301;
t326 = t295 * t302;
t325 = qJD(2) * t299;
t324 = qJD(2) * t302;
t323 = qJD(4) * t298;
t321 = t295 * t324;
t320 = t295 * t323;
t319 = t295 * t325;
t318 = t299 * t331;
t316 = t331 * t301;
t315 = t294 * t318;
t314 = t296 * t317;
t282 = t294 * t299 - t314;
t310 = t282 * t301 + t296 * t329;
t309 = t284 * t298 + t294 * t327;
t307 = -pkin(2) - pkin(5) - pkin(8) - t313;
t283 = t294 * t302 + t296 * t318;
t286 = t298 * t331 + t301 * t326;
t303 = qJD(3) - t305 * t301 + (t298 * t311 + t301 * t322) * qJD(4);
t285 = t296 * t302 - t315;
t281 = -qJD(2) * t315 + t296 * t324;
t280 = t284 * qJD(2);
t279 = t283 * qJD(2);
t278 = -qJD(2) * t314 + t294 * t325;
t272 = qJD(4) * t316 - t301 * t319 - t302 * t320;
t267 = t294 * t329 - t330;
t266 = -t282 * t323 + (qJD(4) * t295 * t296 + t279) * t301;
t264 = qJD(4) * t309 - t281 * t301;
t1 = [0, -t280 * t332 + t281 * t307 + t284 * t308 + t285 * t303, t281, t305 * t309 - t322 * t264 - t311 * (-qJD(4) * t330 - t281 * t298 + t294 * t320) t264 (t264 * t300 + t280 * t297) * r_i_i_C(1) + (-t264 * t297 + t280 * t300) * r_i_i_C(2) + ((-t267 * t297 - t285 * t300) * r_i_i_C(1) + (-t267 * t300 + t285 * t297) * r_i_i_C(2)) * qJD(6); 0, -t278 * t332 + t279 * t307 + t282 * t308 + t283 * t303, t279, -t305 * (-t282 * t298 + t296 * t327) + t322 * t266 + t311 * (qJD(4) * t310 + t279 * t298) -t266 (-t266 * t300 + t278 * t297) * r_i_i_C(1) + (t266 * t297 + t278 * t300) * r_i_i_C(2) + ((-t283 * t300 + t297 * t310) * r_i_i_C(1) + (t283 * t297 + t300 * t310) * r_i_i_C(2)) * qJD(6); 0 ((t332 * qJD(2) - t308) * t302 + (qJD(2) * t307 + t303) * t299) * t295, t319, t305 * (-t298 * t326 + t316) - t322 * t272 - t311 * (qJD(4) * t286 - t298 * t319) t272 (t272 * t300 - t297 * t321) * r_i_i_C(1) + (-t272 * t297 - t300 * t321) * r_i_i_C(2) + ((-t286 * t297 - t300 * t328) * r_i_i_C(1) + (-t286 * t300 + t297 * t328) * r_i_i_C(2)) * qJD(6);];
JaD_transl  = t1;
