% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP3
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
% Datum: 2019-02-26 20:57
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:36
% EndTime: 2019-02-26 20:57:36
% DurationCPUTime: 0.20s
% Computational Cost: add. (318->54), mult. (514->84), div. (0->0), fcn. (429->8), ass. (0->41)
t248 = sin(qJ(3));
t250 = cos(qJ(3));
t280 = pkin(8) + r_i_i_C(1);
t282 = t280 * t250;
t284 = (-pkin(3) * t248 + t282) * qJD(3);
t247 = sin(qJ(4));
t249 = cos(qJ(4));
t277 = r_i_i_C(3) + qJ(5);
t279 = r_i_i_C(2) - pkin(4);
t256 = -t277 * t247 + t279 * t249;
t253 = -pkin(3) + t256;
t252 = t253 * t248 + t282;
t255 = t279 * t247 + t277 * t249;
t281 = t255 * qJD(4) + qJD(5) * t247;
t246 = qJ(1) + pkin(9);
t245 = cos(t246);
t276 = t245 * t247;
t275 = t247 * t250;
t274 = t249 * t250;
t244 = sin(t246);
t273 = qJD(1) * t244;
t272 = qJD(1) * t245;
t271 = qJD(3) * t248;
t270 = qJD(3) * t250;
t269 = qJD(4) * t247;
t268 = qJD(4) * t249;
t266 = t280 * t248;
t263 = t249 * t273;
t262 = t244 * t271;
t261 = t244 * t269;
t260 = t245 * t268;
t259 = t244 * t247 + t245 * t274;
t258 = t244 * t275 + t245 * t249;
t257 = -pkin(3) * t250 - pkin(2) - t266;
t254 = t244 * t268 + t247 * t272;
t251 = -t281 * t248 + (t253 * t250 - t266) * qJD(3);
t233 = t259 * qJD(1) - t249 * t262 - t250 * t261 - t260;
t232 = -t245 * t269 - t247 * t262 + t254 * t250 - t263;
t231 = t250 * t263 + (t249 * t271 + t250 * t269) * t245 - t254;
t230 = t258 * qJD(1) - t250 * t260 + t271 * t276 - t261;
t1 = [-t258 * qJD(5) + t279 * t233 - t277 * t232 - t244 * t284 + (-cos(qJ(1)) * pkin(1) - t244 * pkin(7) + t257 * t245) * qJD(1), 0, t251 * t245 - t252 * t273, t259 * qJD(5) - t279 * t230 - t277 * t231, -t230, 0; -(t244 * t249 - t245 * t275) * qJD(5) + t279 * t231 - t277 * t230 + t245 * t284 + (-sin(qJ(1)) * pkin(1) + t245 * pkin(7) + t257 * t244) * qJD(1), 0, t251 * t244 + t252 * t272 -(-t244 * t274 + t276) * qJD(5) + t277 * t233 + t279 * t232, t232, 0; 0, 0, t252 * qJD(3) + t281 * t250, t255 * t270 + (t256 * qJD(4) + t249 * qJD(5)) * t248, t247 * t270 + t248 * t268, 0;];
JaD_transl  = t1;
