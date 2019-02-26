% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:35
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:35:23
% EndTime: 2019-02-26 21:35:23
% DurationCPUTime: 0.23s
% Computational Cost: add. (308->56), mult. (536->83), div. (0->0), fcn. (446->8), ass. (0->45)
t262 = qJ(2) + pkin(9);
t260 = sin(t262);
t261 = cos(t262);
t304 = pkin(8) + r_i_i_C(2);
t281 = t304 * t261 - sin(qJ(2)) * pkin(2);
t309 = -pkin(3) * t260 + t281;
t264 = sin(qJ(4));
t267 = cos(qJ(4));
t299 = r_i_i_C(3) + qJ(5);
t303 = r_i_i_C(1) + pkin(4);
t276 = -t299 * t264 - t303 * t267;
t272 = -pkin(3) + t276;
t271 = t272 * t260 + t281;
t305 = t303 * t264 - t299 * t267;
t308 = t305 * qJD(4) - qJD(5) * t264;
t306 = -t304 * t260 - cos(qJ(2)) * pkin(2);
t266 = sin(qJ(1));
t298 = t266 * t264;
t297 = t266 * t267;
t269 = cos(qJ(1));
t296 = t269 * t264;
t295 = t269 * t267;
t294 = qJD(1) * t266;
t293 = qJD(1) * t269;
t292 = qJD(2) * t261;
t291 = qJD(2) * t266;
t290 = qJD(2) * t269;
t289 = qJD(4) * t267;
t288 = qJD(4) * t269;
t285 = t260 * t291;
t284 = qJD(4) * t298;
t283 = t260 * t290;
t282 = t267 * t288;
t279 = t261 * t295 + t298;
t278 = t261 * t298 + t295;
t277 = -pkin(3) * t261 - pkin(1) + t306;
t274 = t264 * t288 + t267 * t294;
t273 = t264 * t293 + t266 * t289;
t270 = t308 * t260 + (t272 * t261 + t306) * qJD(2);
t263 = -qJ(3) - pkin(7);
t248 = t279 * qJD(1) - t261 * t284 - t267 * t285 - t282;
t247 = t273 * t261 - t264 * t285 - t274;
t246 = t274 * t261 + t267 * t283 - t273;
t245 = t278 * qJD(1) - t261 * t282 + t264 * t283 - t284;
t1 = [-t278 * qJD(5) + t269 * qJD(3) - t303 * t248 - t299 * t247 - t309 * t291 + (t266 * t263 + t277 * t269) * qJD(1), t270 * t269 - t271 * t294, t293, t279 * qJD(5) + t303 * t245 - t299 * t246, -t245, 0; -(-t261 * t296 + t297) * qJD(5) + t266 * qJD(3) - t303 * t246 - t299 * t245 + t309 * t290 + (-t269 * t263 + t277 * t266) * qJD(1), t270 * t266 + t271 * t293, t294 -(-t261 * t297 + t296) * qJD(5) + t299 * t248 - t303 * t247, t247, 0; 0, t271 * qJD(2) - t308 * t261, 0, -t305 * t292 + (t276 * qJD(4) + t267 * qJD(5)) * t260, t260 * t289 + t264 * t292, 0;];
JaD_transl  = t1;
