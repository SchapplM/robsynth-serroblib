% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:36
% EndTime: 2019-02-26 22:10:37
% DurationCPUTime: 0.30s
% Computational Cost: add. (402->68), mult. (437->101), div. (0->0), fcn. (338->10), ass. (0->56)
t260 = pkin(10) + qJ(5);
t256 = sin(t260);
t257 = cos(t260);
t262 = qJ(2) + qJ(3);
t258 = sin(t262);
t297 = qJD(5) * t258;
t259 = cos(t262);
t261 = qJD(2) + qJD(3);
t304 = t259 * t261;
t319 = t256 * t304 + t257 * t297;
t254 = cos(pkin(10)) * pkin(4) + pkin(3);
t318 = r_i_i_C(1) * t257 + t254;
t253 = t258 * qJD(4);
t265 = sin(qJ(2));
t306 = pkin(2) * qJD(2);
t295 = t265 * t306;
t264 = -pkin(9) - qJ(4);
t307 = r_i_i_C(3) - t264;
t317 = (-t254 * t258 + t307 * t259) * t261 + (pkin(4) * sin(pkin(10)) + pkin(8) + pkin(7)) * qJD(1) + t253 - t295;
t289 = t256 * t297;
t316 = r_i_i_C(1) * t289 + t319 * r_i_i_C(2) + qJD(4) * t259;
t268 = cos(qJ(1));
t282 = qJD(5) * t259 - qJD(1);
t315 = t268 * t282;
t281 = qJD(1) * t259 - qJD(5);
t266 = sin(qJ(1));
t302 = t261 * t266;
t293 = t258 * t302;
t313 = t281 * t268 - t293;
t311 = pkin(2) * t265;
t309 = r_i_i_C(2) * t256;
t308 = r_i_i_C(3) * t259;
t303 = t259 * t264;
t301 = t261 * t268;
t300 = qJD(1) * t266;
t299 = qJD(1) * t268;
t296 = t258 * t309;
t292 = t258 * t301;
t291 = t258 * t300;
t290 = t258 * t299;
t280 = t318 * t261;
t279 = t282 * t266;
t277 = t264 * t293 + t316 * t266 + t290 * t309 + t299 * t308;
t276 = -t258 * t318 - t303;
t275 = t264 * t292 + t316 * t268 + t318 * t291 + t300 * t303;
t267 = cos(qJ(2));
t274 = qJD(1) * (-t267 * pkin(2) - t254 * t259 - t307 * t258 - pkin(1));
t273 = (-r_i_i_C(3) * t258 - t259 * t318) * t261;
t272 = t281 * t266 + t292;
t271 = -t267 * t306 + t273;
t270 = t261 * t296 + r_i_i_C(3) * t304 + t253 - t258 * t280 + (-t261 * t264 + (-r_i_i_C(1) * t256 - r_i_i_C(2) * t257) * qJD(5)) * t259;
t235 = t256 * t279 - t313 * t257;
t234 = t313 * t256 + t257 * t279;
t233 = t256 * t315 + t272 * t257;
t232 = t272 * t256 - t257 * t315;
t1 = [t235 * r_i_i_C(1) + t234 * r_i_i_C(2) - t317 * t266 + t268 * t274 (-t296 - t308 + t311) * t300 + t271 * t268 + t275 (-r_i_i_C(3) * t301 - t300 * t309) * t258 + (-r_i_i_C(3) * t300 - t268 * t280) * t259 + t275, t259 * t301 - t291, t232 * r_i_i_C(1) + t233 * r_i_i_C(2), 0; -t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t266 * t274 + t317 * t268, t271 * t266 + (t276 - t311) * t299 + t277, t266 * t273 + t276 * t299 + t277, t259 * t302 + t290, -t234 * r_i_i_C(1) + t235 * r_i_i_C(2), 0; 0, t270 - t295, t270, t261 * t258 (-t257 * t304 + t289) * r_i_i_C(2) - t319 * r_i_i_C(1), 0;];
JaD_transl  = t1;
