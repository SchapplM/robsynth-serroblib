% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:03
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:03:14
% EndTime: 2019-02-26 21:03:14
% DurationCPUTime: 0.29s
% Computational Cost: add. (478->64), mult. (467->88), div. (0->0), fcn. (354->9), ass. (0->56)
t262 = pkin(10) + qJ(3);
t260 = qJ(4) + t262;
t257 = cos(t260);
t264 = sin(qJ(6));
t266 = cos(qJ(6));
t312 = r_i_i_C(1) * t264 + r_i_i_C(2) * t266 + qJ(5);
t316 = t257 * t312;
t263 = qJD(3) + qJD(4);
t303 = t257 * t263;
t252 = qJ(5) * t303;
t256 = sin(t260);
t294 = pkin(4) + pkin(9) + r_i_i_C(3);
t283 = t294 * t263;
t258 = sin(t262);
t304 = pkin(3) * qJD(3);
t293 = t258 * t304;
t315 = (qJD(5) - t283) * t256 + (pkin(5) + pkin(8) + pkin(7) + qJ(2)) * qJD(1) + t252 - t293;
t295 = qJD(6) * t266;
t286 = t257 * t295;
t314 = r_i_i_C(1) * t286 + qJD(5) * t257;
t281 = qJD(6) * t256 + qJD(1);
t311 = t266 * t281;
t310 = t281 * t264;
t308 = pkin(3) * t258;
t302 = t263 * t264;
t301 = t263 * t266;
t267 = cos(qJ(1));
t300 = t263 * t267;
t265 = sin(qJ(1));
t299 = qJD(1) * t265;
t298 = qJD(1) * t267;
t296 = qJD(6) * t264;
t292 = t257 * t301;
t291 = t265 * t303;
t290 = t257 * t300;
t288 = t256 * t299;
t287 = t257 * t296;
t285 = t294 * t256;
t284 = t294 * t257;
t280 = -qJD(1) * t256 - qJD(6);
t279 = t314 * t265 + t298 * t316;
t278 = t314 * t267 + t294 * t288;
t277 = t280 * t267;
t274 = t312 * t256;
t273 = -r_i_i_C(2) * t296 - t283;
t272 = t280 * t265 + t290;
t259 = cos(t262);
t271 = qJD(2) + (-qJ(5) * t256 - pkin(3) * t259 - cos(pkin(10)) * pkin(2) - pkin(1) - t284) * qJD(1);
t270 = t257 * r_i_i_C(1) * t302 + r_i_i_C(2) * t292 + t252 + (r_i_i_C(1) * t295 + qJD(5) + t273) * t256;
t269 = -r_i_i_C(2) * t287 + (-t284 - t274) * t263;
t268 = -t259 * t304 + t269;
t239 = t272 * t264 + t267 * t311;
t238 = t272 * t266 - t267 * t310;
t237 = -t265 * t311 + (t277 - t291) * t264;
t236 = t266 * t277 + (-t292 + t310) * t265;
t1 = [t237 * r_i_i_C(1) + t236 * r_i_i_C(2) - t315 * t265 + t271 * t267, t298 (t308 - t316) * t299 + t268 * t267 + t278, -t274 * t300 + (t273 * t267 - t299 * t312) * t257 + t278, -t288 + t290, t238 * r_i_i_C(1) - t239 * r_i_i_C(2); t239 * r_i_i_C(1) + t238 * r_i_i_C(2) + t271 * t265 + t315 * t267, t299 (-t285 - t308) * t298 + t268 * t265 + t279, t269 * t265 - t285 * t298 + t279, t256 * t298 + t291, -t236 * r_i_i_C(1) + t237 * r_i_i_C(2); 0, 0, t270 - t293, t270, t263 * t256 (-t256 * t302 + t286) * r_i_i_C(2) + (t256 * t301 + t287) * r_i_i_C(1);];
JaD_transl  = t1;
