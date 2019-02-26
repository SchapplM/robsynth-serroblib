% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_transl = S6RPRRPP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:57:27
% EndTime: 2019-02-26 20:57:27
% DurationCPUTime: 0.29s
% Computational Cost: add. (432->62), mult. (688->93), div. (0->0), fcn. (584->8), ass. (0->43)
t250 = sin(qJ(3));
t252 = cos(qJ(3));
t270 = pkin(5) + pkin(8) + r_i_i_C(1);
t283 = t270 * t252;
t287 = (-pkin(3) * t250 + t283) * qJD(3);
t249 = sin(qJ(4));
t251 = cos(qJ(4));
t269 = pkin(4) + r_i_i_C(3) + qJ(6);
t279 = r_i_i_C(2) + qJ(5);
t281 = t269 * t249 - t279 * t251;
t286 = -t270 * qJD(3) + t281 * qJD(4) - qJD(5) * t249 - qJD(6) * t251;
t257 = -t279 * t249 - t269 * t251;
t255 = -pkin(3) + t257;
t284 = t255 * t250 + t283;
t278 = t249 * t252;
t277 = t251 * t252;
t248 = qJ(1) + pkin(9);
t246 = sin(t248);
t276 = qJD(1) * t246;
t247 = cos(t248);
t275 = qJD(1) * t247;
t274 = qJD(3) * t250;
t273 = qJD(3) * t252;
t272 = qJD(4) * t249;
t271 = qJD(4) * t251;
t268 = t251 * t276;
t267 = t249 * t274;
t266 = t251 * t274;
t265 = t246 * t272;
t264 = t247 * t271;
t260 = t246 * t249 + t247 * t277;
t232 = t246 * t278 + t247 * t251;
t259 = -pkin(3) * t252 - t270 * t250 - pkin(2);
t258 = t246 * t271 + t249 * t275;
t254 = qJD(3) * t255;
t253 = t286 * t250 + t252 * t254;
t234 = -t246 * t251 + t247 * t278;
t233 = t246 * t277 - t247 * t249;
t231 = t260 * qJD(1) - t246 * t266 - t252 * t265 - t264;
t230 = -t246 * t267 - t247 * t272 + t258 * t252 - t268;
t229 = t252 * t268 + (t252 * t272 + t266) * t247 - t258;
t228 = t232 * qJD(1) + t247 * t267 - t252 * t264 - t265;
t1 = [-t232 * qJD(5) - t233 * qJD(6) - t279 * t230 - t269 * t231 - t246 * t287 + (-cos(qJ(1)) * pkin(1) - pkin(7) * t246 + t259 * t247) * qJD(1), 0, t253 * t247 - t284 * t276, qJD(5) * t260 - t234 * qJD(6) + t269 * t228 - t279 * t229, -t228, -t229; t234 * qJD(5) + t260 * qJD(6) - t279 * t228 - t269 * t229 + t247 * t287 + (-sin(qJ(1)) * pkin(1) + pkin(7) * t247 + t259 * t246) * qJD(1), 0, t253 * t246 + t284 * t275, t233 * qJD(5) - t232 * qJD(6) - t269 * t230 + t279 * t231, t230, t231; 0, 0, t250 * t254 - t286 * t252, -t281 * t273 + (t257 * qJD(4) + qJD(5) * t251 - qJD(6) * t249) * t250, t249 * t273 + t250 * t271, -t250 * t272 + t251 * t273;];
JaD_transl  = t1;
