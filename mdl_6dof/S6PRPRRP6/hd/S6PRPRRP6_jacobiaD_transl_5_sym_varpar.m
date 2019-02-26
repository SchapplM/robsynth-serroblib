% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRP6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP6_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:06
% EndTime: 2019-02-26 19:53:07
% DurationCPUTime: 0.31s
% Computational Cost: add. (236->64), mult. (762->121), div. (0->0), fcn. (760->10), ass. (0->49)
t292 = sin(pkin(10));
t294 = cos(pkin(10));
t298 = sin(qJ(2));
t295 = cos(pkin(6));
t301 = cos(qJ(2));
t321 = t295 * t301;
t284 = t292 * t321 + t294 * t298;
t300 = cos(qJ(4));
t293 = sin(pkin(6));
t297 = sin(qJ(4));
t326 = t293 * t297;
t330 = -t284 * t300 + t292 * t326;
t296 = sin(qJ(5));
t299 = cos(qJ(5));
t313 = r_i_i_C(1) * t299 - r_i_i_C(2) * t296;
t305 = t313 * qJD(5);
t311 = pkin(4) + t313;
t328 = pkin(9) + r_i_i_C(3);
t329 = t311 * t297 - t328 * t300 + qJ(3);
t325 = t293 * t298;
t324 = t293 * t300;
t323 = t293 * t301;
t322 = t295 * t298;
t320 = qJD(2) * t298;
t319 = qJD(2) * t301;
t317 = t292 * t320;
t316 = t293 * t319;
t315 = t294 * t319;
t314 = t293 * t320;
t312 = -r_i_i_C(1) * t296 - r_i_i_C(2) * t299;
t310 = -pkin(2) - pkin(8) + t312;
t282 = t292 * t298 - t294 * t321;
t309 = -t282 * t297 + t294 * t324;
t308 = t282 * t300 + t294 * t326;
t271 = t284 * t297 + t292 * t324;
t283 = t292 * t301 + t294 * t322;
t307 = t295 * t297 + t300 * t323;
t306 = -t295 * t300 + t297 * t323;
t304 = qJD(5) * t312;
t302 = qJD(3) + t297 * t304 + (t328 * t297 + t311 * t300) * qJD(4);
t285 = -t292 * t322 + t294 * t301;
t281 = -t295 * t317 + t315;
t280 = t284 * qJD(2);
t279 = t283 * qJD(2);
t278 = -t295 * t315 + t317;
t274 = t307 * qJD(4) - t297 * t314;
t268 = t308 * qJD(4) + t279 * t297;
t266 = t330 * qJD(4) - t281 * t297;
t1 = [0, -t280 * t329 + t310 * t281 - t284 * t305 + t302 * t285, t281, -t328 * t266 - t330 * t304 + t311 * (-t271 * qJD(4) + t281 * t300) (t266 * t296 - t280 * t299) * r_i_i_C(1) + (t266 * t299 + t280 * t296) * r_i_i_C(2) + ((-t271 * t299 - t285 * t296) * r_i_i_C(1) + (t271 * t296 - t285 * t299) * r_i_i_C(2)) * qJD(5), 0; 0, -t278 * t329 + t310 * t279 - t282 * t305 + t302 * t283, t279, t328 * t268 + t308 * t304 + t311 * (t309 * qJD(4) + t279 * t300) (-t268 * t296 - t278 * t299) * r_i_i_C(1) + (-t268 * t299 + t278 * t296) * r_i_i_C(2) + ((-t283 * t296 + t299 * t309) * r_i_i_C(1) + (-t283 * t299 - t296 * t309) * r_i_i_C(2)) * qJD(5), 0; 0 ((t329 * qJD(2) + t305) * t301 + (t310 * qJD(2) + t302) * t298) * t293, t314, -t328 * t274 - t307 * t304 + t311 * (t306 * qJD(4) + t300 * t314) (t274 * t296 + t299 * t316) * r_i_i_C(1) + (t274 * t299 - t296 * t316) * r_i_i_C(2) + ((-t296 * t325 + t299 * t306) * r_i_i_C(1) + (-t296 * t306 - t299 * t325) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
