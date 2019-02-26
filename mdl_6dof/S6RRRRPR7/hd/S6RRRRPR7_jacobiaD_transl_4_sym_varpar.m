% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRPR7
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR7_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR7_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR7_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR7_jacobiaD_transl_4_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:34:07
% EndTime: 2019-02-26 22:34:08
% DurationCPUTime: 0.29s
% Computational Cost: add. (332->69), mult. (663->114), div. (0->0), fcn. (618->10), ass. (0->55)
t280 = cos(pkin(6));
t282 = sin(qJ(2));
t283 = sin(qJ(1));
t314 = t283 * t282;
t302 = t280 * t314;
t285 = cos(qJ(2));
t286 = cos(qJ(1));
t311 = t286 * t285;
t266 = -qJD(1) * t302 - qJD(2) * t314 + (qJD(2) * t280 + qJD(1)) * t311;
t277 = qJD(3) + qJD(4);
t279 = sin(pkin(6));
t315 = t279 * t286;
t322 = t277 * t315 - t266;
t321 = r_i_i_C(3) + pkin(10) + pkin(9);
t320 = pkin(3) * qJD(3);
t312 = t286 * t282;
t313 = t283 * t285;
t268 = t280 * t312 + t313;
t319 = t268 * t277;
t278 = qJ(3) + qJ(4);
t275 = sin(t278);
t318 = t275 * t277;
t317 = t279 * t283;
t284 = cos(qJ(3));
t316 = t279 * t284;
t276 = cos(t278);
t292 = t302 - t311;
t306 = qJD(1) * t286;
t300 = t279 * t306;
t290 = t277 * t292 + t300;
t293 = t280 * t313 + t312;
t264 = t268 * qJD(1) + t293 * qJD(2);
t296 = t277 * t317 - t264;
t259 = -t296 * t275 + t290 * t276;
t260 = t290 * t275 + t296 * t276;
t310 = t259 * r_i_i_C(1) - t260 * r_i_i_C(2);
t307 = qJD(1) * t283;
t301 = t279 * t307;
t291 = t301 - t319;
t298 = t322 * t276;
t309 = (t322 * t275 + t291 * t276) * r_i_i_C(1) + (-t291 * t275 + t298) * r_i_i_C(2);
t299 = qJD(2) * t279 * t285;
t289 = -t277 * t280 - t299;
t304 = t277 * t279 * t282;
t308 = (t289 * t275 - t276 * t304) * r_i_i_C(1) + (t275 * t304 + t289 * t276) * r_i_i_C(2);
t281 = sin(qJ(3));
t305 = t281 * t320;
t297 = -r_i_i_C(1) * t275 - r_i_i_C(2) * t276;
t274 = t284 * pkin(3) + pkin(2);
t295 = t276 * r_i_i_C(1) - t275 * r_i_i_C(2) + t274;
t294 = t280 * t311 - t314;
t288 = t297 * t277 - t305;
t265 = t293 * qJD(1) + t268 * qJD(2);
t263 = -t294 * qJD(1) + t292 * qJD(2);
t1 = [(t268 * t318 + t298) * r_i_i_C(1) + (t266 * t275 + t276 * t319) * r_i_i_C(2) - t266 * t274 + t268 * t305 - pkin(1) * t306 - t321 * t265 + ((-r_i_i_C(2) * t318 + t284 * t320) * t286 + (-pkin(3) * t281 - pkin(8) + t297) * t307) * t279, t295 * t263 - t321 * t264 - t288 * t293 (t284 * t300 + t264 * t281 + (-t281 * t317 + t284 * t292) * qJD(3)) * pkin(3) + t310, t310, 0, 0; t260 * r_i_i_C(1) + t259 * r_i_i_C(2) - t264 * t274 - t321 * t263 + (-pkin(1) * t283 + pkin(8) * t315) * qJD(1) + (t281 * t300 + (t281 * t292 + t283 * t316) * qJD(3)) * pkin(3), -t295 * t265 + t321 * t266 + t288 * t294 (t284 * t301 - t266 * t281 + (-t268 * t284 + t281 * t315) * qJD(3)) * pkin(3) + t309, t309, 0, 0; 0 (t288 * t285 + (-t295 * t282 + t321 * t285) * qJD(2)) * t279 (-t281 * t299 + (-t280 * t281 - t282 * t316) * qJD(3)) * pkin(3) + t308, t308, 0, 0;];
JaD_transl  = t1;
