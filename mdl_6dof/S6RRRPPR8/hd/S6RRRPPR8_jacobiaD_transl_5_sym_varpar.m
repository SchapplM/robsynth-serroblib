% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR8
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:07
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPPR8_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPPR8_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPPR8_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPPR8_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:07:40
% EndTime: 2019-02-26 22:07:41
% DurationCPUTime: 0.31s
% Computational Cost: add. (360->66), mult. (1060->106), div. (0->0), fcn. (1028->8), ass. (0->46)
t265 = sin(qJ(3));
t268 = cos(qJ(3));
t284 = pkin(3) + pkin(4) - r_i_i_C(2);
t297 = r_i_i_C(1) + qJ(4);
t301 = (t284 * t265 - t297 * t268) * qJD(3) - qJD(4) * t265;
t264 = cos(pkin(6));
t267 = sin(qJ(1));
t266 = sin(qJ(2));
t290 = t267 * t266;
t280 = t264 * t290;
t287 = qJD(2) * t266;
t269 = cos(qJ(2));
t270 = cos(qJ(1));
t288 = t269 * t270;
t250 = -qJD(1) * t280 - t267 * t287 + (qJD(2) * t264 + qJD(1)) * t288;
t289 = t267 * t269;
t291 = t266 * t270;
t255 = t264 * t291 + t289;
t263 = sin(pkin(6));
t293 = t263 * t270;
t277 = t255 * t265 + t268 * t293;
t295 = t263 * t267;
t282 = t265 * t295;
t300 = -qJD(1) * t282 + t277 * qJD(3) - t250 * t268;
t298 = t297 * t265 + t284 * t268 + pkin(2);
t274 = t280 - t288;
t296 = t274 * t265;
t294 = t263 * t268;
t292 = t265 * t270;
t286 = qJD(3) * t268;
t283 = pkin(9) - r_i_i_C(3) - qJ(5);
t281 = t263 * t292;
t279 = qJD(1) * t294;
t278 = qJD(2) * t263 * t269;
t276 = -t268 * t274 + t282;
t275 = t264 * t265 + t266 * t294;
t254 = t264 * t288 - t290;
t256 = t264 * t289 + t291;
t243 = -qJD(3) * t281 + t250 * t265 + t255 * t286 - t267 * t279;
t251 = t275 * qJD(3) + t265 * t278;
t249 = t256 * qJD(1) + t255 * qJD(2);
t248 = t255 * qJD(1) + t256 * qJD(2);
t247 = -t254 * qJD(1) + t274 * qJD(2);
t242 = -t248 * t268 + qJD(3) * t296 + (qJD(1) * t292 + t267 * t286) * t263;
t241 = t276 * qJD(3) - t248 * t265 - t270 * t279;
t1 = [-t254 * qJD(5) - t277 * qJD(4) - t250 * pkin(2) - t297 * t243 + (-pkin(1) * t270 - pkin(8) * t295) * qJD(1) - t283 * t249 + t284 * t300, qJD(5) * t274 + t298 * t247 - t283 * t248 + t301 * t256, t276 * qJD(4) - t284 * t241 + t297 * t242, t241, t247, 0; -t256 * qJD(5) - (t267 * t294 + t296) * qJD(4) - t248 * pkin(2) + t297 * t241 + (-pkin(1) * t267 + pkin(8) * t293) * qJD(1) - t283 * t247 + t284 * t242, -qJD(5) * t255 - t249 * t298 + t283 * t250 - t254 * t301 -(-t255 * t268 + t281) * qJD(4) - t297 * t300 - t284 * t243, t243, -t249, 0; 0 (-qJD(5) * t266 - t301 * t269 + (-t266 * t298 + t283 * t269) * qJD(2)) * t263, t275 * qJD(4) + t297 * (t268 * t278 + (-t263 * t265 * t266 + t264 * t268) * qJD(3)) - t284 * t251, t251, -t263 * t287, 0;];
JaD_transl  = t1;
