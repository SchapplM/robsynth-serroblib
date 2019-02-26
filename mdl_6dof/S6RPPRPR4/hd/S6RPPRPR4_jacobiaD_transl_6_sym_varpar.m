% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:27
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR4_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR4_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:27:25
% EndTime: 2019-02-26 20:27:25
% DurationCPUTime: 0.25s
% Computational Cost: add. (265->55), mult. (532->83), div. (0->0), fcn. (510->10), ass. (0->41)
t248 = qJ(4) + pkin(10);
t246 = sin(t248);
t247 = cos(t248);
t278 = pkin(8) + r_i_i_C(3);
t259 = -t278 * t247 + sin(qJ(4)) * pkin(4);
t250 = sin(qJ(6));
t252 = cos(qJ(6));
t260 = t252 * r_i_i_C(1) - t250 * r_i_i_C(2) + pkin(5);
t280 = t260 * t246 + t259;
t285 = t280 * qJD(4);
t266 = t278 * t246;
t284 = -t260 * t247 - t266;
t272 = sin(pkin(9));
t273 = cos(pkin(9));
t276 = sin(qJ(1));
t277 = cos(qJ(1));
t240 = t277 * t272 - t276 * t273;
t282 = qJD(1) * t276;
t281 = qJD(1) * t277;
t279 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t274 = cos(qJ(4)) * pkin(4);
t271 = qJD(4) * t246;
t270 = qJD(4) * t247;
t269 = qJD(6) * t247;
t268 = qJD(6) * t250;
t267 = qJD(6) * t252;
t263 = t250 * r_i_i_C(1) + t252 * r_i_i_C(2);
t239 = -t276 * t272 - t277 * t273;
t237 = t239 * qJD(1);
t262 = t239 * t269 + t237;
t238 = t240 * qJD(1);
t261 = t240 * t269 + t238;
t258 = qJD(6) * t263;
t257 = qJD(6) * t239 + t237 * t247 - t240 * t271;
t256 = qJD(6) * t240 + t238 * t247 + t239 * t271;
t254 = -t246 * t258 + (t274 - t284) * qJD(4);
t249 = -qJ(5) - pkin(7);
t245 = pkin(3) + t274;
t236 = t250 * t262 + t256 * t252;
t235 = -t256 * t250 + t252 * t262;
t1 = [(-t238 * t250 + t239 * t267) * r_i_i_C(1) + (-t238 * t252 - t239 * t268) * r_i_i_C(2) + t238 * t249 + t239 * qJD(5) - qJ(2) * t282 - (-t245 + t284) * t237 + (-t247 * t258 - t285) * t240 + t279 * t277, t281, 0, -t238 * t280 + t254 * t239, t237, t235 * r_i_i_C(1) - t236 * r_i_i_C(2); t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t240 * qJD(5) - t237 * t249 + qJ(2) * t281 + (pkin(5) * t247 + t245 + t266) * t238 + (pkin(5) * t246 + t259) * t239 * qJD(4) + t279 * t276, t282, 0, t237 * t280 + t254 * t240, t238 (r_i_i_C(1) * t261 + t257 * r_i_i_C(2)) * t252 + (t257 * r_i_i_C(1) - r_i_i_C(2) * t261) * t250; 0, 0, 0, t263 * t269 + t285, 0 (-t246 * t268 + t252 * t270) * r_i_i_C(2) + (t246 * t267 + t250 * t270) * r_i_i_C(1);];
JaD_transl  = t1;
