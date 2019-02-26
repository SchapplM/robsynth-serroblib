% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:05
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:05:10
% EndTime: 2019-02-26 22:05:10
% DurationCPUTime: 0.34s
% Computational Cost: add. (354->67), mult. (610->99), div. (0->0), fcn. (506->8), ass. (0->51)
t263 = cos(qJ(3));
t301 = t263 * pkin(3);
t255 = pkin(2) + t301;
t260 = sin(qJ(3));
t261 = sin(qJ(2));
t264 = cos(qJ(2));
t258 = qJ(3) + pkin(10);
t256 = sin(t258);
t300 = r_i_i_C(2) + qJ(4) + pkin(8);
t273 = t300 * qJD(2) + t256 * qJD(5);
t298 = pkin(3) * qJD(3);
t302 = pkin(3) * t260;
t311 = (-t260 * t298 + t273) * t264 + (pkin(7) + t302) * qJD(1) - (qJD(2) * t255 - qJD(4)) * t261;
t257 = cos(t258);
t299 = r_i_i_C(3) + qJ(5);
t303 = -r_i_i_C(1) - pkin(4);
t269 = t303 * t256 + t299 * t257 - t302;
t309 = qJD(3) * t269 + t273;
t272 = -t299 * t256 + t303 * t257;
t270 = -t255 + t272;
t308 = t270 * t261 + t300 * t264;
t262 = sin(qJ(1));
t297 = t262 * t264;
t265 = cos(qJ(1));
t296 = t265 * t257;
t295 = qJD(1) * t262;
t294 = qJD(1) * t264;
t293 = qJD(1) * t265;
t292 = qJD(2) * t261;
t291 = qJD(2) * t264;
t290 = qJD(2) * t265;
t289 = qJD(3) * t262;
t288 = qJD(3) * t265;
t286 = t257 * qJD(5);
t283 = t262 * t292;
t282 = t261 * t290;
t281 = t256 * t289;
t280 = t256 * t288;
t279 = t257 * t288;
t278 = -qJD(3) + t294;
t276 = (-qJD(3) * t264 + qJD(1)) * t263;
t275 = t262 * t256 + t264 * t296;
t271 = t256 * t293 + t257 * t289;
t268 = t270 * qJD(2) + qJD(4);
t267 = t263 * t298 - t286 + (-t255 * t264 - t300 * t261 - pkin(1)) * qJD(1);
t266 = -t309 * t261 + t268 * t264;
t244 = t275 * qJD(1) - t257 * t283 - t264 * t281 - t279;
t243 = -t256 * t283 - t257 * t295 + t271 * t264 - t280;
t242 = t264 * t280 + (t262 * t294 + t282) * t257 - t271;
t241 = t256 * t282 - t264 * t279 - t281 + (t256 * t297 + t296) * qJD(1);
t1 = [-t299 * t243 + t303 * t244 - t311 * t262 + t267 * t265, t266 * t265 - t308 * t295, t275 * qJD(5) - t299 * t242 - t303 * t241 + (t265 * t276 + (t262 * t278 + t282) * t260) * pkin(3), -t261 * t295 + t264 * t290, -t241, 0; -t299 * t241 + t303 * t242 + t267 * t262 + t311 * t265, t266 * t262 + t308 * t293 -(t265 * t256 - t257 * t297) * qJD(5) + t299 * t244 + t303 * t243 + (t262 * t276 + (-t265 * t278 + t283) * t260) * pkin(3), t261 * t293 + t262 * t291, t243, 0; 0, t268 * t261 + t309 * t264, t269 * t291 + (t286 + (t272 - t301) * qJD(3)) * t261, t292, t261 * qJD(3) * t257 + t256 * t291, 0;];
JaD_transl  = t1;
