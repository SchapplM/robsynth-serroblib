% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:49
% EndTime: 2019-02-26 20:58:50
% DurationCPUTime: 0.22s
% Computational Cost: add. (301->58), mult. (516->89), div. (0->0), fcn. (433->7), ass. (0->44)
t261 = sin(qJ(4));
t263 = cos(qJ(4));
t292 = r_i_i_C(3) + qJ(5);
t295 = pkin(4) + r_i_i_C(1);
t296 = t295 * t261 - t292 * t263;
t298 = -t296 * qJD(4) + qJD(5) * t261;
t259 = pkin(9) + qJ(3);
t257 = sin(t259);
t258 = cos(t259);
t294 = pkin(8) + r_i_i_C(2);
t279 = t294 * t258;
t297 = -pkin(3) * t257 + t279;
t271 = -t292 * t261 - t295 * t263;
t267 = -pkin(3) + t271;
t262 = sin(qJ(1));
t291 = t262 * t261;
t290 = t262 * t263;
t264 = cos(qJ(1));
t289 = t264 * t261;
t288 = t264 * t263;
t287 = qJD(1) * t262;
t286 = qJD(1) * t264;
t285 = qJD(3) * t258;
t284 = qJD(3) * t262;
t283 = qJD(3) * t264;
t282 = qJD(4) * t263;
t281 = qJD(4) * t264;
t278 = t257 * t284;
t277 = qJD(4) * t291;
t276 = t257 * t283;
t275 = t263 * t281;
t274 = t258 * t288 + t291;
t273 = t258 * t291 + t288;
t272 = -pkin(3) * t258 - t294 * t257 - cos(pkin(9)) * pkin(2) - pkin(1);
t269 = t261 * t281 + t263 * t287;
t268 = t261 * t286 + t262 * t282;
t266 = qJD(3) * t267;
t265 = -t294 * qJD(3) - t298;
t260 = -pkin(7) - qJ(2);
t245 = t274 * qJD(1) - t258 * t277 - t263 * t278 - t275;
t244 = t268 * t258 - t261 * t278 - t269;
t243 = t269 * t258 + t263 * t276 - t268;
t242 = t273 * qJD(1) - t258 * t275 + t261 * t276 - t277;
t1 = [-t273 * qJD(5) + t264 * qJD(2) - t295 * t245 - t292 * t244 - t297 * t284 + (t262 * t260 + t272 * t264) * qJD(1), t286 (t264 * t266 - t294 * t287) * t258 + (t265 * t264 - t267 * t287) * t257, t274 * qJD(5) + t295 * t242 - t292 * t243, -t242, 0; -(-t258 * t289 + t290) * qJD(5) + t262 * qJD(2) - t295 * t243 - t292 * t242 + t297 * t283 + (-t264 * t260 + t272 * t262) * qJD(1), t287 (t262 * t266 + t294 * t286) * t258 + (t265 * t262 + t267 * t286) * t257 -(-t258 * t290 + t289) * qJD(5) + t292 * t245 - t295 * t244, t244, 0; 0, 0, t298 * t258 + (t267 * t257 + t279) * qJD(3), -t296 * t285 + (t271 * qJD(4) + t263 * qJD(5)) * t257, t257 * t282 + t261 * t285, 0;];
JaD_transl  = t1;
