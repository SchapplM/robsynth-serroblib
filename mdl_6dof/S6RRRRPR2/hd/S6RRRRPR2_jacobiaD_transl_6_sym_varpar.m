% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:31
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRPR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRPR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRPR2_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:31:17
% EndTime: 2019-02-26 22:31:17
% DurationCPUTime: 0.33s
% Computational Cost: add. (759->82), mult. (577->110), div. (0->0), fcn. (437->12), ass. (0->66)
t279 = pkin(11) + qJ(6);
t273 = sin(t279);
t274 = cos(t279);
t282 = qJ(2) + qJ(3);
t278 = qJ(4) + t282;
t271 = sin(t278);
t321 = qJD(6) * t271;
t272 = cos(t278);
t280 = qJD(2) + qJD(3);
t275 = qJD(4) + t280;
t328 = t272 * t275;
t344 = t273 * t328 + t274 * t321;
t270 = cos(pkin(11)) * pkin(5) + pkin(4);
t343 = r_i_i_C(1) * t274 + t270;
t269 = t271 * qJD(5);
t285 = sin(qJ(2));
t330 = pkin(2) * qJD(2);
t317 = t285 * t330;
t276 = sin(t282);
t335 = pkin(3) * t280;
t320 = t276 * t335;
t284 = -pkin(10) - qJ(5);
t331 = r_i_i_C(3) - t284;
t342 = (-t270 * t271 + t331 * t272) * t275 + (pkin(5) * sin(pkin(11)) + pkin(7) + pkin(8) + pkin(9)) * qJD(1) + t269 - t317 - t320;
t311 = t273 * t321;
t341 = r_i_i_C(1) * t311 + t344 * r_i_i_C(2) + qJD(5) * t272;
t288 = cos(qJ(1));
t303 = qJD(6) * t272 - qJD(1);
t340 = t288 * t303;
t277 = cos(t282);
t294 = (-r_i_i_C(3) * t271 - t272 * t343) * t275;
t291 = -t277 * t335 + t294;
t302 = qJD(1) * t272 - qJD(6);
t286 = sin(qJ(1));
t326 = t275 * t286;
t316 = t271 * t326;
t338 = t302 * t288 - t316;
t336 = pkin(3) * t276;
t333 = r_i_i_C(2) * t273;
t332 = r_i_i_C(3) * t272;
t327 = t272 * t284;
t325 = t275 * t288;
t324 = qJD(1) * t286;
t323 = qJD(1) * t288;
t318 = t271 * t333;
t315 = t271 * t325;
t313 = t271 * t324;
t312 = t271 * t323;
t301 = t343 * t275;
t300 = t303 * t286;
t299 = -t318 - t332;
t298 = t284 * t316 + t341 * t286 + t312 * t333 + t323 * t332;
t297 = -t271 * t343 - t327;
t296 = t284 * t315 + t341 * t288 + t343 * t313 + t324 * t327;
t287 = cos(qJ(2));
t295 = qJD(1) * (-t287 * pkin(2) - pkin(3) * t277 - t270 * t272 - t331 * t271 - pkin(1));
t293 = t302 * t286 + t315;
t292 = -t287 * t330 + t291;
t290 = t275 * t318 + r_i_i_C(3) * t328 + t269 - t271 * t301 + (-t275 * t284 + (-r_i_i_C(1) * t273 - r_i_i_C(2) * t274) * qJD(6)) * t272;
t289 = t290 - t320;
t266 = -t285 * pkin(2) - t336;
t247 = t273 * t300 - t338 * t274;
t246 = t338 * t273 + t274 * t300;
t245 = t273 * t340 + t293 * t274;
t244 = t293 * t273 - t274 * t340;
t1 = [t247 * r_i_i_C(1) + t246 * r_i_i_C(2) - t342 * t286 + t288 * t295 (-t266 + t299) * t324 + t292 * t288 + t296 (t299 + t336) * t324 + t291 * t288 + t296 (-r_i_i_C(3) * t325 - t324 * t333) * t271 + (-r_i_i_C(3) * t324 - t288 * t301) * t272 + t296, t272 * t325 - t313, t244 * r_i_i_C(1) + t245 * r_i_i_C(2); -t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t286 * t295 + t342 * t288, t292 * t286 + (t266 + t297) * t323 + t298, t291 * t286 + (t297 - t336) * t323 + t298, t286 * t294 + t297 * t323 + t298, t272 * t326 + t312, -t246 * r_i_i_C(1) + t247 * r_i_i_C(2); 0, t289 - t317, t289, t290, t275 * t271 (-t274 * t328 + t311) * r_i_i_C(2) - t344 * r_i_i_C(1);];
JaD_transl  = t1;
