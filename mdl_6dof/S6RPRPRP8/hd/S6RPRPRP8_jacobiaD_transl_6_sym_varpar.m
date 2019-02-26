% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:47:37
% EndTime: 2019-02-26 20:47:37
% DurationCPUTime: 0.29s
% Computational Cost: add. (313->55), mult. (548->80), div. (0->0), fcn. (454->8), ass. (0->43)
t248 = qJ(3) + pkin(9);
t245 = sin(t248);
t246 = cos(t248);
t250 = sin(qJ(5));
t253 = cos(qJ(5));
t281 = r_i_i_C(3) + qJ(6);
t286 = -r_i_i_C(1) - pkin(5);
t262 = -t281 * t250 + t286 * t253;
t260 = -pkin(4) + t262;
t287 = pkin(8) + r_i_i_C(2);
t266 = t287 * t246 - sin(qJ(3)) * pkin(3);
t274 = qJD(6) * t250;
t288 = t286 * t250 + t281 * t253;
t290 = t288 * qJD(5) + t274;
t295 = t290 * t246 + (t260 * t245 + t266) * qJD(3);
t293 = -pkin(4) * t245 - qJ(2) + t266;
t264 = t287 * t245 + cos(qJ(3)) * pkin(3);
t291 = t260 * t246 - t264;
t282 = -pkin(1) - qJ(4) - pkin(7);
t252 = sin(qJ(1));
t280 = t252 * t250;
t279 = t252 * t253;
t255 = cos(qJ(1));
t278 = t255 * t250;
t277 = qJD(1) * t252;
t247 = qJD(1) * t255;
t276 = qJD(3) * t245;
t275 = qJD(3) * t246;
t273 = t253 * qJD(6);
t272 = t255 * t245 * t253;
t271 = t255 * t275;
t270 = qJD(5) * t245 + qJD(1);
t269 = qJD(1) * t245 + qJD(5);
t268 = -qJD(4) + t273;
t267 = t269 * t255;
t263 = t245 * t279 + t278;
t257 = t245 * t274 + qJD(2) + (pkin(4) * t246 + t264) * qJD(3);
t256 = qJD(1) * t291;
t240 = t253 * t267 + (-t270 * t250 + t253 * t275) * t252;
t239 = t270 * t279 + (t252 * t275 + t267) * t250;
t238 = -t253 * t271 + (t245 * t278 + t279) * qJD(5) + t263 * qJD(1);
t237 = -qJD(5) * t272 - t253 * t247 - t250 * t271 + t269 * t280;
t1 = [t268 * t252 + t286 * t238 - t281 * t237 + t257 * t255 + (t293 * t252 + t282 * t255) * qJD(1), t247, t252 * t295 - t255 * t256, -t277, t263 * qJD(6) + t286 * t239 + t281 * t240, t239; -t268 * t255 - t286 * t240 + t281 * t239 + t257 * t252 + (t282 * t252 - t293 * t255) * qJD(1), t277, -t252 * t256 - t255 * t295, t247 -(t272 - t280) * qJD(6) + t281 * t238 + t286 * t237, t237; 0, 0, t291 * qJD(3) - t245 * t290, 0, -t288 * t276 + (t262 * qJD(5) + t273) * t246, t246 * qJD(5) * t253 - t250 * t276;];
JaD_transl  = t1;
