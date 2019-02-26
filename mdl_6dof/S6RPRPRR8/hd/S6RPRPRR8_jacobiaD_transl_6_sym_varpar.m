% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:52
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:52:47
% EndTime: 2019-02-26 20:52:48
% DurationCPUTime: 0.28s
% Computational Cost: add. (379->56), mult. (458->79), div. (0->0), fcn. (359->10), ass. (0->50)
t253 = cos(qJ(5));
t240 = t253 * pkin(5) + pkin(4);
t247 = qJ(3) + pkin(10);
t241 = sin(t247);
t242 = cos(t247);
t285 = r_i_i_C(3) + pkin(9) + pkin(8);
t263 = t285 * t242 - sin(qJ(3)) * pkin(3);
t276 = qJD(5) * t253;
t306 = -pkin(5) * t276 + (-t240 * t241 - qJ(2) + t263) * qJD(1) - qJD(4);
t246 = qJD(5) + qJD(6);
t250 = sin(qJ(5));
t290 = pkin(5) * t250;
t275 = qJD(5) * t290;
t248 = qJ(5) + qJ(6);
t244 = sin(t248);
t245 = cos(t248);
t301 = r_i_i_C(1) * t244 + r_i_i_C(2) * t245;
t259 = t301 * t246 + t275;
t300 = -r_i_i_C(1) * t245 + r_i_i_C(2) * t244;
t264 = t240 - t300;
t305 = -t259 * t242 + (-t264 * t241 + t263) * qJD(3);
t261 = t285 * t241 + cos(qJ(3)) * pkin(3);
t302 = t264 * t242 + t261;
t299 = t253 * (qJD(5) * t241 + qJD(1));
t252 = sin(qJ(1));
t280 = qJD(1) * t241;
t271 = t246 + t280;
t255 = cos(qJ(1));
t277 = qJD(3) * t242;
t273 = t255 * t277;
t294 = t271 * t252 - t273;
t274 = t252 * t277;
t293 = t271 * t255 + t274;
t272 = -t241 * t246 - qJD(1);
t265 = t272 * t255;
t233 = t294 * t244 + t245 * t265;
t234 = t244 * t265 - t294 * t245;
t282 = -t233 * r_i_i_C(1) + t234 * r_i_i_C(2);
t266 = t272 * t252;
t235 = -t293 * t244 + t245 * t266;
t236 = t244 * t266 + t293 * t245;
t281 = t235 * r_i_i_C(1) - t236 * r_i_i_C(2);
t279 = qJD(1) * t252;
t278 = qJD(3) * t241;
t269 = -qJD(5) - t280;
t260 = t300 * t242 * t246 + t301 * t278;
t258 = qJD(1) * t302;
t257 = -t241 * t275 + qJD(2) + (-pkin(1) - qJ(4) - pkin(7) - t290) * qJD(1) + (t240 * t242 + t261) * qJD(3);
t243 = qJD(1) * t255;
t1 = [t234 * r_i_i_C(1) + t233 * r_i_i_C(2) + t306 * t252 + t257 * t255, t243, t305 * t252 + t255 * t258, -t279 (-t252 * t299 + (t269 * t255 - t274) * t250) * pkin(5) + t281, t281; t236 * r_i_i_C(1) + t235 * r_i_i_C(2) + t257 * t252 - t306 * t255, t279, t252 * t258 - t305 * t255, t243 (t255 * t299 + (t269 * t252 + t273) * t250) * pkin(5) + t282, t282; 0, 0, -t302 * qJD(3) + t259 * t241, 0 (-t242 * t276 + t250 * t278) * pkin(5) + t260, t260;];
JaD_transl  = t1;
