% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:53
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR1_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR1_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:53:44
% EndTime: 2019-02-26 19:53:44
% DurationCPUTime: 0.24s
% Computational Cost: add. (275->52), mult. (634->101), div. (0->0), fcn. (640->12), ass. (0->48)
t283 = sin(pkin(12));
t286 = cos(pkin(12));
t292 = cos(qJ(2));
t305 = qJD(2) * t292;
t290 = sin(qJ(2));
t306 = qJD(2) * t290;
t317 = t283 * t306 - t286 * t305;
t316 = -r_i_i_C(3) - pkin(9) - pkin(8);
t315 = pkin(2) * qJD(2);
t282 = qJ(4) + qJ(5);
t279 = sin(t282);
t281 = qJD(4) + qJD(5);
t314 = t279 * t281;
t280 = cos(t282);
t313 = t280 * t281;
t285 = sin(pkin(6));
t312 = t281 * t285;
t289 = sin(qJ(4));
t311 = t285 * t289;
t288 = cos(pkin(6));
t310 = t288 * t290;
t298 = t292 * t283 + t290 * t286;
t269 = t298 * t288;
t284 = sin(pkin(11));
t287 = cos(pkin(11));
t297 = t290 * t283 - t292 * t286;
t259 = t287 * t269 - t284 * t297;
t264 = t317 * t288;
t271 = -t283 * t305 - t286 * t306;
t299 = t287 * t264 - t284 * t271;
t301 = t287 * t312 + t299;
t309 = (-t259 * t313 + t301 * t279) * r_i_i_C(1) + (t259 * t314 + t301 * t280) * r_i_i_C(2);
t261 = -t284 * t269 - t287 * t297;
t256 = t284 * t264 + t287 * t271;
t300 = -t284 * t312 - t256;
t308 = (-t261 * t313 + t300 * t279) * r_i_i_C(1) + (t261 * t314 + t300 * t280) * r_i_i_C(2);
t267 = t298 * t285;
t262 = t317 * t285;
t302 = -t281 * t288 + t262;
t307 = (-t267 * t313 + t302 * t279) * r_i_i_C(1) + (t267 * t314 + t302 * t280) * r_i_i_C(2);
t291 = cos(qJ(4));
t296 = t291 * pkin(4) + r_i_i_C(1) * t280 - r_i_i_C(2) * t279 + pkin(3);
t295 = qJD(2) * t298;
t294 = -pkin(4) * qJD(4) * t289 + (-r_i_i_C(1) * t279 - r_i_i_C(2) * t280) * t281;
t270 = t297 * qJD(2);
t268 = t297 * t288;
t265 = t288 * t295;
t1 = [0, -t316 * t256 + (t284 * t310 - t287 * t292) * t315 + t294 * (t284 * t268 - t287 * t298) + t296 * (t284 * t265 + t287 * t270) 0 (-t256 * t289 + (-t261 * t291 - t284 * t311) * qJD(4)) * pkin(4) + t308, t308, 0; 0, t316 * t299 + (-t284 * t292 - t287 * t310) * t315 + t294 * (-t287 * t268 - t284 * t298) + t296 * (-t287 * t265 + t284 * t270) 0 (t299 * t289 + (-t259 * t291 + t287 * t311) * qJD(4)) * pkin(4) + t309, t309, 0; 0, t316 * t262 + (-pkin(2) * t306 - t294 * t297 - t295 * t296) * t285, 0 (t262 * t289 + (-t267 * t291 - t288 * t289) * qJD(4)) * pkin(4) + t307, t307, 0;];
JaD_transl  = t1;
