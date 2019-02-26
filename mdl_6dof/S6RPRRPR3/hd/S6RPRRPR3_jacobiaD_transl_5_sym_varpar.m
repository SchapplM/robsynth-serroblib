% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR3
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
% Datum: 2019-02-26 21:02
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:02:03
% EndTime: 2019-02-26 21:02:03
% DurationCPUTime: 0.21s
% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
t250 = sin(qJ(3));
t252 = cos(qJ(3));
t282 = pkin(8) + r_i_i_C(2);
t284 = t282 * t252;
t286 = (-pkin(3) * t250 + t284) * qJD(3);
t249 = sin(qJ(4));
t251 = cos(qJ(4));
t279 = r_i_i_C(3) + qJ(5);
t281 = -r_i_i_C(1) - pkin(4);
t258 = -t279 * t249 + t281 * t251;
t255 = -pkin(3) + t258;
t254 = t255 * t250 + t284;
t257 = t281 * t249 + t279 * t251;
t283 = t257 * qJD(4) + qJD(5) * t249;
t248 = qJ(1) + pkin(10);
t247 = cos(t248);
t278 = t247 * t249;
t277 = t249 * t252;
t276 = t251 * t252;
t246 = sin(t248);
t275 = qJD(1) * t246;
t274 = qJD(1) * t247;
t273 = qJD(3) * t250;
t272 = qJD(3) * t252;
t271 = qJD(4) * t249;
t270 = qJD(4) * t251;
t268 = t282 * t250;
t265 = t251 * t275;
t264 = t246 * t273;
t263 = t246 * t271;
t262 = t247 * t270;
t261 = t246 * t249 + t247 * t276;
t260 = t246 * t277 + t247 * t251;
t259 = -pkin(3) * t252 - pkin(2) - t268;
t256 = t246 * t270 + t249 * t274;
t253 = -t283 * t250 + (t255 * t252 - t268) * qJD(3);
t235 = t261 * qJD(1) - t251 * t264 - t252 * t263 - t262;
t234 = -t247 * t271 - t249 * t264 + t256 * t252 - t265;
t233 = t252 * t265 + (t251 * t273 + t252 * t271) * t247 - t256;
t232 = t260 * qJD(1) - t252 * t262 + t273 * t278 - t263;
t1 = [-t260 * qJD(5) + t281 * t235 - t279 * t234 - t246 * t286 + (-cos(qJ(1)) * pkin(1) - t246 * pkin(7) + t259 * t247) * qJD(1), 0, t253 * t247 - t254 * t275, t261 * qJD(5) - t281 * t232 - t279 * t233, -t232, 0; -(t246 * t251 - t247 * t277) * qJD(5) + t281 * t233 - t279 * t232 + t247 * t286 + (-sin(qJ(1)) * pkin(1) + t247 * pkin(7) + t259 * t246) * qJD(1), 0, t253 * t246 + t254 * t274 -(-t246 * t276 + t278) * qJD(5) + t279 * t235 + t281 * t234, t234, 0; 0, 0, t254 * qJD(3) + t283 * t252, t257 * t272 + (t258 * qJD(4) + t251 * qJD(5)) * t250, t249 * t272 + t250 * t270, 0;];
JaD_transl  = t1;
