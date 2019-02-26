% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:04
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:04:33
% EndTime: 2019-02-26 22:04:34
% DurationCPUTime: 0.32s
% Computational Cost: add. (410->65), mult. (545->89), div. (0->0), fcn. (404->8), ass. (0->55)
t264 = qJ(2) + qJ(3);
t262 = cos(t264);
t268 = cos(qJ(6));
t311 = pkin(5) + qJ(4);
t319 = r_i_i_C(1) * t268 + t311;
t320 = t262 * t319;
t261 = sin(t264);
t259 = t261 * qJD(4);
t263 = qJD(2) + qJD(3);
t297 = pkin(3) + pkin(4) + pkin(9) + r_i_i_C(3);
t286 = t297 * t261;
t266 = sin(qJ(2));
t310 = pkin(2) * qJD(2);
t298 = t266 * t310;
t318 = (-t262 * t311 + t286) * t263 + (qJ(5) - pkin(8) - pkin(7)) * qJD(1) - t259 + t298;
t265 = sin(qJ(6));
t308 = t263 * t261;
t296 = t265 * t308;
t316 = r_i_i_C(2) * t296 + qJD(4) * t262;
t313 = pkin(2) * t266;
t309 = t262 * t263;
t307 = t263 * t268;
t270 = cos(qJ(1));
t306 = t263 * t270;
t305 = t268 * t270;
t267 = sin(qJ(1));
t303 = qJD(1) * t267;
t302 = qJD(1) * t270;
t300 = qJD(6) * t262;
t299 = r_i_i_C(2) * t262 * t265;
t295 = t262 * t307;
t294 = t267 * t309;
t293 = t262 * t306;
t291 = t261 * t303;
t289 = qJD(6) * t261 + qJD(1);
t288 = qJD(1) * t261 + qJD(6);
t285 = t297 * t262;
t284 = t316 * t267 + t302 * t320;
t283 = t289 * t265;
t282 = t319 * t261;
t280 = (-r_i_i_C(1) * t265 - r_i_i_C(2) * t268) * qJD(6);
t279 = t316 * t270 + t297 * t291 + t299 * t303;
t278 = -t286 - t299;
t277 = t267 * t288 - t293;
t276 = -t263 * t297 + t280;
t269 = cos(qJ(2));
t275 = -qJD(5) + (-t269 * pkin(2) - t261 * t311 - pkin(1) - t285) * qJD(1);
t274 = t262 * t280 + (-t282 - t285) * t263;
t273 = r_i_i_C(1) * t295 + t261 * t276 - t263 * t299 + t311 * t309 + t259;
t272 = -t269 * t310 + t274;
t243 = -t288 * t305 + (t283 - t295) * t267;
t242 = t289 * t268 * t267 + (t270 * t288 + t294) * t265;
t241 = t268 * t277 + t270 * t283;
t240 = t265 * t277 - t289 * t305;
t1 = [t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t318 * t267 + t275 * t270 (t313 - t320) * t303 + t272 * t270 + t279, -t282 * t306 + (t270 * t276 - t303 * t319) * t262 + t279, -t291 + t293, -t302, t240 * r_i_i_C(1) + t241 * r_i_i_C(2); -t241 * r_i_i_C(1) + t240 * r_i_i_C(2) + t275 * t267 - t318 * t270 (t278 - t313) * t302 + t272 * t267 + t284, t267 * t274 + t278 * t302 + t284, t261 * t302 + t294, -t303, -t242 * r_i_i_C(1) + t243 * r_i_i_C(2); 0, t273 - t298, t273, t308, 0 (-t261 * t307 - t265 * t300) * r_i_i_C(2) + (t268 * t300 - t296) * r_i_i_C(1);];
JaD_transl  = t1;
