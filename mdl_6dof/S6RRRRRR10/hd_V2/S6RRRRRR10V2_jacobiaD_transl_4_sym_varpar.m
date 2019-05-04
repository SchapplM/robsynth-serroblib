% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRRRR10V2
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
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-11 14:56
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR10V2_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR10V2_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S6RRRRRR10V2_jacobiaD_transl_4_sym_varpar: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-11 14:56:37
% EndTime: 2019-04-11 14:56:38
% DurationCPUTime: 0.30s
% Computational Cost: add. (266->55), mult. (390->92), div. (0->0), fcn. (295->8), ass. (0->53)
t257 = qJ(2) + qJ(3);
t254 = sin(t257);
t261 = cos(qJ(4));
t308 = r_i_i_C(1) * t261 + pkin(3);
t311 = t254 * t308;
t258 = sin(qJ(4));
t291 = qJD(4) * t261;
t255 = cos(t257);
t256 = qJD(2) + qJD(3);
t299 = t255 * t256;
t310 = t254 * t291 + t258 * t299;
t305 = pkin(5) + r_i_i_C(3);
t286 = t305 * t255;
t259 = sin(qJ(2));
t300 = pkin(2) * qJD(2);
t288 = t259 * t300;
t303 = pkin(3) * t254;
t309 = -t288 + (t286 - t303) * t256;
t292 = qJD(4) * t258;
t279 = t254 * t292;
t306 = r_i_i_C(1) * t279 + t310 * r_i_i_C(2);
t304 = pkin(2) * t259;
t301 = r_i_i_C(2) * t258;
t260 = sin(qJ(1));
t298 = t256 * t260;
t297 = t256 * t261;
t263 = cos(qJ(1));
t296 = t256 * t263;
t295 = t261 * t263;
t294 = qJD(1) * t260;
t293 = qJD(1) * t263;
t290 = t254 * t301;
t289 = qJD(1) * t301;
t287 = t305 * t254;
t285 = t305 * t260;
t284 = t254 * t297;
t274 = qJD(4) * t255 - qJD(1);
t273 = qJD(1) * t255 - qJD(4);
t272 = t308 * t255;
t271 = t308 * t263;
t270 = t306 * t263 + t294 * t311;
t269 = t274 * t258;
t268 = t263 * t254 * t289 + t306 * t260 + t293 * t286;
t262 = cos(qJ(2));
t267 = qJD(1) * (-t262 * pkin(2) - pkin(3) * t255 - pkin(1) - t287);
t266 = t254 * t296 + t273 * t260;
t265 = -t262 * t300 + (-t272 - t287) * t256;
t264 = -t255 * r_i_i_C(2) * t291 + (-t255 * t292 - t284) * r_i_i_C(1) + t305 * t299 + (-t303 + t290) * t256;
t238 = -t273 * t295 + (t269 + t284) * t260;
t237 = t274 * t261 * t260 + (-t254 * t298 + t273 * t263) * t258;
t236 = t266 * t261 + t263 * t269;
t235 = t266 * t258 - t274 * t295;
t1 = [t238 * r_i_i_C(1) + t237 * r_i_i_C(2) - t309 * t260 + t263 * t267 (-t286 - t290 + t304) * t294 + t265 * t263 + t270 (-t260 * t289 - t305 * t296) * t254 + (-qJD(1) * t285 - t256 * t271) * t255 + t270, t235 * r_i_i_C(1) + t236 * r_i_i_C(2), 0, 0; -t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t260 * t267 + t309 * t263 (-t304 - t311) * t293 + t265 * t260 + t268, -t272 * t298 + (-qJD(1) * t271 - t256 * t285) * t254 + t268, -t237 * r_i_i_C(1) + t238 * r_i_i_C(2), 0, 0; 0, t264 - t288, t264 (-t255 * t297 + t279) * r_i_i_C(2) - t310 * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
