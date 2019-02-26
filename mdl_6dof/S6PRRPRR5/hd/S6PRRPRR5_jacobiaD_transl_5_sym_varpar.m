% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRRPRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRRPRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRRPRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRRPRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:06:38
% EndTime: 2019-02-26 20:06:39
% DurationCPUTime: 0.37s
% Computational Cost: add. (336->69), mult. (819->125), div. (0->0), fcn. (824->12), ass. (0->51)
t304 = sin(qJ(3));
t306 = cos(qJ(3));
t299 = pkin(12) + qJ(5);
t297 = sin(t299);
t298 = cos(t299);
t319 = t297 * r_i_i_C(1) + t298 * r_i_i_C(2);
t314 = qJD(5) * t319;
t320 = t298 * r_i_i_C(1) - t297 * r_i_i_C(2);
t317 = cos(pkin(12)) * pkin(4) + pkin(3) + t320;
t336 = r_i_i_C(3) + pkin(9) + qJ(4);
t308 = (t317 * t304 - t336 * t306) * qJD(3) - t304 * qJD(4) + t306 * t314;
t301 = sin(pkin(11));
t305 = sin(qJ(2));
t307 = cos(qJ(2));
t334 = cos(pkin(11));
t335 = cos(pkin(6));
t318 = t335 * t334;
t287 = t301 * t307 + t305 * t318;
t302 = sin(pkin(6));
t324 = t302 * t334;
t277 = t287 * t306 - t304 * t324;
t325 = t301 * t335;
t289 = -t305 * t325 + t334 * t307;
t332 = t302 * t304;
t331 = t302 * t306;
t330 = t302 * t307;
t329 = qJD(2) * t305;
t327 = t302 * t329;
t326 = qJD(2) * t330;
t316 = t307 * t318;
t315 = -t289 * t304 + t301 * t331;
t279 = t289 * t306 + t301 * t332;
t313 = t320 * qJD(5);
t312 = -t287 * t304 - t306 * t324;
t311 = -t305 * t332 + t335 * t306;
t291 = t335 * t304 + t305 * t331;
t310 = sin(pkin(12)) * pkin(4) + pkin(8) + t319;
t288 = t334 * t305 + t307 * t325;
t309 = -t336 * t304 - t317 * t306 - pkin(2);
t286 = t301 * t305 - t316;
t285 = t289 * qJD(2);
t284 = t288 * qJD(2);
t283 = t287 * qJD(2);
t282 = -qJD(2) * t316 + t301 * t329;
t281 = t311 * qJD(3) + t306 * t326;
t280 = t291 * qJD(3) + t304 * t326;
t275 = t315 * qJD(3) - t284 * t306;
t274 = t279 * qJD(3) - t284 * t304;
t273 = t312 * qJD(3) - t282 * t306;
t272 = t277 * qJD(3) - t282 * t304;
t1 = [0, -t310 * t284 + t309 * t285 + t308 * t288 + t289 * t313, t279 * qJD(4) - t317 * t274 + t336 * t275 - t315 * t314, t274 (-t275 * t297 + t285 * t298) * r_i_i_C(1) + (-t275 * t298 - t285 * t297) * r_i_i_C(2) + ((-t279 * t298 - t288 * t297) * r_i_i_C(1) + (t279 * t297 - t288 * t298) * r_i_i_C(2)) * qJD(5), 0; 0, -t310 * t282 + t309 * t283 + t308 * t286 + t287 * t313, t277 * qJD(4) - t317 * t272 + t336 * t273 - t312 * t314, t272 (-t273 * t297 + t283 * t298) * r_i_i_C(1) + (-t273 * t298 - t283 * t297) * r_i_i_C(2) + ((-t277 * t298 - t286 * t297) * r_i_i_C(1) + (t277 * t297 - t286 * t298) * r_i_i_C(2)) * qJD(5), 0; 0 ((t309 * qJD(2) + t313) * t305 + (t310 * qJD(2) - t308) * t307) * t302, t291 * qJD(4) - t317 * t280 + t336 * t281 - t311 * t314, t280 (-t281 * t297 + t298 * t327) * r_i_i_C(1) + (-t281 * t298 - t297 * t327) * r_i_i_C(2) + ((-t291 * t298 + t297 * t330) * r_i_i_C(1) + (t291 * t297 + t298 * t330) * r_i_i_C(2)) * qJD(5), 0;];
JaD_transl  = t1;
