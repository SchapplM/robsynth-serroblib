% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:26
% EndTime: 2019-02-26 20:59:26
% DurationCPUTime: 0.33s
% Computational Cost: add. (358->62), mult. (618->92), div. (0->0), fcn. (512->8), ass. (0->47)
t245 = sin(qJ(3));
t248 = cos(qJ(3));
t247 = cos(qJ(4));
t281 = t247 * pkin(4);
t239 = pkin(3) + t281;
t242 = qJ(4) + pkin(9);
t240 = sin(t242);
t241 = cos(t242);
t279 = r_i_i_C(3) + qJ(6);
t283 = -r_i_i_C(1) - pkin(5);
t257 = -t279 * t240 + t283 * t241;
t286 = -t239 + t257;
t253 = t286 * qJD(3) + qJD(5);
t280 = r_i_i_C(2) + qJ(5) + pkin(8);
t258 = t280 * qJD(3) + t240 * qJD(6);
t244 = sin(qJ(4));
t282 = pkin(4) * t244;
t284 = t283 * t240 + t279 * t241 - t282;
t288 = t284 * qJD(4) + t258;
t290 = t253 * t245 + t288 * t248;
t268 = t241 * qJD(6);
t278 = pkin(4) * qJD(4);
t287 = (-t239 * t245 + t280 * t248 - qJ(2)) * qJD(1) - t247 * t278 + t268;
t246 = sin(qJ(1));
t249 = cos(qJ(1));
t262 = qJD(1) * t245 + qJD(4);
t271 = qJD(3) * t248;
t255 = t246 * t271 + t262 * t249;
t277 = t246 * t240;
t276 = t246 * t241;
t275 = t249 * t240;
t274 = qJD(1) * t246;
t273 = qJD(1) * t249;
t272 = qJD(3) * t245;
t270 = qJD(3) * t249;
t267 = t249 * t245 * t241;
t265 = t248 * t270;
t263 = qJD(4) * t245 + qJD(1);
t261 = t263 * t246;
t259 = t245 * t276 + t275;
t252 = qJD(1) * (t280 * t245 - t248 * t286);
t250 = qJD(2) + (qJD(3) * t239 - qJD(5)) * t248 + (-pkin(1) - pkin(7) - t282) * qJD(1) + (-t244 * t278 + t258) * t245;
t234 = -t240 * t261 + t255 * t241;
t233 = t255 * t240 + t263 * t276;
t232 = -t241 * t265 + (t245 * t275 + t276) * qJD(4) + t259 * qJD(1);
t231 = -qJD(4) * t267 - t240 * t265 - t241 * t273 + t262 * t277;
t1 = [-t279 * t231 + t283 * t232 + t287 * t246 + t250 * t249, t273, t290 * t246 + t249 * t252, t259 * qJD(6) + t279 * t234 + t283 * t233 + (-t255 * t244 - t247 * t261) * pkin(4), t246 * t272 - t248 * t273, t233; t279 * t233 - t283 * t234 + t250 * t246 - t287 * t249, t274, t246 * t252 - t290 * t249 -(t267 - t277) * qJD(6) + t279 * t232 + t283 * t231 + (t263 * t249 * t247 + (-t262 * t246 + t265) * t244) * pkin(4), -t245 * t270 - t248 * t274, t231; 0, 0, -t245 * t288 + t253 * t248, -t284 * t272 + (t268 + (t257 - t281) * qJD(4)) * t248, t271, t248 * qJD(4) * t241 - t240 * t272;];
JaD_transl  = t1;
