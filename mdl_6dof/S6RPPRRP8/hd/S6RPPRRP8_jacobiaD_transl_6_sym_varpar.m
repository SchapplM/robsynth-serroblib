% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP8_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:34:13
% EndTime: 2019-02-26 20:34:13
% DurationCPUTime: 0.21s
% Computational Cost: add. (306->57), mult. (528->87), div. (0->0), fcn. (441->7), ass. (0->44)
t248 = sin(qJ(5));
t250 = cos(qJ(5));
t278 = r_i_i_C(3) + qJ(6);
t281 = pkin(5) + r_i_i_C(1);
t258 = -t278 * t248 - t281 * t250;
t256 = -pkin(4) + t258;
t245 = pkin(9) + qJ(4);
t242 = sin(t245);
t280 = pkin(8) + r_i_i_C(2);
t284 = t280 * t242;
t243 = cos(t245);
t283 = t280 * t243 - pkin(3) * sin(pkin(9)) - pkin(4) * t242 - qJ(2);
t257 = t281 * t248 - t278 * t250;
t269 = qJD(6) * t248;
t253 = t257 * qJD(5) - t269;
t282 = -t280 * qJD(4) + t253;
t279 = -pkin(1) - pkin(7) - qJ(3);
t277 = t243 * t250;
t249 = sin(qJ(1));
t276 = t249 * t248;
t275 = t249 * t250;
t251 = cos(qJ(1));
t274 = t251 * t248;
t273 = qJD(1) * t249;
t244 = qJD(1) * t251;
t272 = qJD(4) * t242;
t271 = qJD(4) * t249;
t270 = qJD(4) * t251;
t268 = t250 * qJD(6);
t267 = t251 * t242 * t250;
t266 = qJD(1) * t280;
t265 = t243 * t270;
t264 = qJD(5) * t242 + qJD(1);
t263 = qJD(1) * t242 + qJD(5);
t262 = -qJD(3) + t268;
t261 = t263 * t251;
t259 = t242 * t275 + t274;
t254 = qJD(1) * t256;
t252 = t242 * t269 + qJD(2) + (pkin(4) * t243 + t284) * qJD(4);
t237 = t250 * t261 + (qJD(4) * t277 - t264 * t248) * t249;
t236 = t264 * t275 + (t243 * t271 + t261) * t248;
t235 = -t250 * t265 + (t242 * t274 + t275) * qJD(5) + t259 * qJD(1);
t234 = -qJD(5) * t267 - t250 * t244 - t248 * t265 + t263 * t276;
t1 = [t262 * t249 - t281 * t235 - t278 * t234 + t252 * t251 + (t283 * t249 + t279 * t251) * qJD(1), t244, -t273 (t251 * t266 + t256 * t271) * t242 + (-t282 * t249 - t251 * t254) * t243, t259 * qJD(6) - t281 * t236 + t278 * t237, t236; -t262 * t251 + t281 * t237 + t278 * t236 + t252 * t249 + (t279 * t249 - t283 * t251) * qJD(1), t273, t244 (t249 * t266 - t256 * t270) * t242 + (-t249 * t254 + t282 * t251) * t243 -(t267 - t276) * qJD(6) + t278 * t235 - t281 * t234, t234; 0, 0, 0, t253 * t242 + (t256 * t243 - t284) * qJD(4), t257 * t272 + (t258 * qJD(5) + t268) * t243, qJD(5) * t277 - t248 * t272;];
JaD_transl  = t1;
