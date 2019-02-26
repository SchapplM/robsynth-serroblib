% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:47
% EndTime: 2019-02-26 20:23:48
% DurationCPUTime: 0.24s
% Computational Cost: add. (252->53), mult. (498->82), div. (0->0), fcn. (479->9), ass. (0->40)
t245 = pkin(10) + qJ(5);
t243 = sin(t245);
t247 = sin(qJ(6));
t248 = cos(qJ(6));
t254 = t248 * r_i_i_C(1) - t247 * r_i_i_C(2) + pkin(5);
t244 = cos(t245);
t271 = pkin(8) + r_i_i_C(3);
t261 = t271 * t244;
t250 = -t254 * t243 + t261;
t278 = qJD(5) * t250;
t260 = t271 * t243;
t277 = t254 * t244 + t260;
t267 = sin(pkin(9));
t268 = cos(pkin(9));
t269 = sin(qJ(1));
t270 = cos(qJ(1));
t237 = t270 * t267 - t269 * t268;
t274 = qJD(1) * t269;
t273 = qJD(1) * t270;
t272 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t266 = qJD(5) * t243;
t265 = qJD(5) * t244;
t264 = qJD(6) * t244;
t263 = qJD(6) * t247;
t262 = qJD(6) * t248;
t257 = t247 * r_i_i_C(1) + t248 * r_i_i_C(2);
t236 = -t267 * t269 - t268 * t270;
t234 = t236 * qJD(1);
t256 = t236 * t264 + t234;
t235 = t237 * qJD(1);
t255 = t237 * t264 + t235;
t253 = qJD(6) * t257;
t252 = qJD(6) * t236 + t234 * t244 - t237 * t266;
t251 = qJD(6) * t237 + t235 * t244 + t236 * t266;
t249 = t277 * qJD(5) - t243 * t253;
t246 = -pkin(7) - qJ(4);
t242 = cos(pkin(10)) * pkin(4) + pkin(3);
t233 = t247 * t256 + t248 * t251;
t232 = -t247 * t251 + t248 * t256;
t1 = [(-t235 * t247 + t236 * t262) * r_i_i_C(1) + (-t235 * t248 - t236 * t263) * r_i_i_C(2) + t235 * t246 + t236 * qJD(4) - qJ(2) * t274 - (-t242 - t277) * t234 + (-t244 * t253 + t278) * t237 + t272 * t270, t273, 0, t234, t235 * t250 + t236 * t249, t232 * r_i_i_C(1) - t233 * r_i_i_C(2); t233 * r_i_i_C(1) + t232 * r_i_i_C(2) + t237 * qJD(4) - t234 * t246 + (pkin(5) * t243 - t261) * t236 * qJD(5) + qJ(2) * t273 + (pkin(5) * t244 + t242 + t260) * t235 + t272 * t269, t274, 0, t235, -t234 * t250 + t237 * t249 (r_i_i_C(1) * t255 + r_i_i_C(2) * t252) * t248 + (r_i_i_C(1) * t252 - r_i_i_C(2) * t255) * t247; 0, 0, 0, 0, t257 * t264 - t278 (-t243 * t263 + t248 * t265) * r_i_i_C(2) + (t243 * t262 + t247 * t265) * r_i_i_C(1);];
JaD_transl  = t1;
