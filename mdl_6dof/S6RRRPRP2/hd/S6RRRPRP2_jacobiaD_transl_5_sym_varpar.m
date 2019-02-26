% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPRP2
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

function JaD_transl = S6RRRPRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:10:05
% EndTime: 2019-02-26 22:10:06
% DurationCPUTime: 0.31s
% Computational Cost: add. (412->63), mult. (426->88), div. (0->0), fcn. (319->10), ass. (0->56)
t271 = qJ(2) + qJ(3);
t266 = pkin(10) + t271;
t264 = sin(t266);
t275 = cos(qJ(5));
t322 = r_i_i_C(1) * t275 + pkin(4);
t289 = t322 * t264;
t272 = sin(qJ(5));
t306 = qJD(5) * t264;
t265 = cos(t266);
t270 = qJD(2) + qJD(3);
t311 = t265 * t270;
t325 = t272 * t311 + t275 * t306;
t319 = pkin(9) + r_i_i_C(3);
t301 = t319 * t265;
t304 = r_i_i_C(2) * t264 * t272;
t324 = t301 + t304;
t267 = sin(t271);
t273 = sin(qJ(2));
t313 = pkin(2) * qJD(2);
t303 = t273 * t313;
t316 = pkin(3) * t270;
t323 = -t267 * t316 + (-pkin(4) * t264 + t301) * t270 - t303;
t296 = t272 * t306;
t320 = r_i_i_C(1) * t296 + t325 * r_i_i_C(2);
t318 = pkin(3) * t267;
t268 = cos(t271);
t317 = pkin(3) * t268;
t312 = t264 * t270;
t310 = t270 * t275;
t277 = cos(qJ(1));
t309 = t275 * t277;
t274 = sin(qJ(1));
t308 = qJD(1) * t274;
t307 = qJD(1) * t277;
t305 = qJD(5) * t265;
t302 = t319 * t264;
t291 = -qJD(1) + t305;
t290 = qJD(1) * t265 - qJD(5);
t288 = t320 * t277 + t308 * t289;
t287 = t291 * t272;
t286 = t320 * t274 + t324 * t307;
t276 = cos(qJ(2));
t284 = -pkin(2) * t276 - pkin(4) * t265 - pkin(1) - t302 - t317;
t283 = -t289 - t318;
t282 = -t265 * t322 - t302;
t281 = t290 * t274 + t277 * t312;
t280 = -t268 * t316 + t282 * t270 - t276 * t313;
t279 = t270 * (t282 - t317);
t278 = (-r_i_i_C(1) * t272 - r_i_i_C(2) * t275) * t305 + t319 * t311 + (t304 + t283) * t270;
t269 = -qJ(4) - pkin(8) - pkin(7);
t263 = -pkin(2) * t273 - t318;
t245 = -t290 * t309 + (t264 * t310 + t287) * t274;
t244 = t291 * t275 * t274 + (-t274 * t312 + t290 * t277) * t272;
t243 = t281 * t275 + t277 * t287;
t242 = t281 * t272 - t291 * t309;
t1 = [t245 * r_i_i_C(1) + t244 * r_i_i_C(2) + t277 * qJD(4) - t323 * t274 + (t269 * t274 + t284 * t277) * qJD(1) (-t263 - t324) * t308 + t280 * t277 + t288 (-t324 + t318) * t308 + t277 * t279 + t288, t307, r_i_i_C(1) * t242 + r_i_i_C(2) * t243, 0; -t243 * r_i_i_C(1) + t242 * r_i_i_C(2) + t274 * qJD(4) + t323 * t277 + (-t269 * t277 + t284 * t274) * qJD(1) (t263 - t289) * t307 + t280 * t274 + t286, t274 * t279 + t283 * t307 + t286, t308, -r_i_i_C(1) * t244 + r_i_i_C(2) * t245, 0; 0, t278 - t303, t278, 0 (-t265 * t310 + t296) * r_i_i_C(2) - t325 * r_i_i_C(1), 0;];
JaD_transl  = t1;
