% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR13
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR13_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR13_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR13_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRR13_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:01:00
% EndTime: 2019-02-26 22:01:00
% DurationCPUTime: 0.19s
% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
t275 = pkin(3) + pkin(8);
t243 = cos(pkin(6));
t246 = sin(qJ(1));
t248 = cos(qJ(2));
t267 = t246 * t248;
t245 = sin(qJ(2));
t249 = cos(qJ(1));
t268 = t245 * t249;
t235 = t243 * t267 + t268;
t251 = t243 * t268 + t267;
t231 = t235 * qJD(1) + t251 * qJD(2);
t244 = sin(qJ(4));
t274 = t231 * t244;
t247 = cos(qJ(4));
t273 = t231 * t247;
t242 = sin(pkin(6));
t272 = t242 * t246;
t271 = t242 * t248;
t270 = t242 * t249;
t269 = t245 * t246;
t266 = t248 * t249;
t265 = qJD(1) * t246;
t264 = qJD(1) * t249;
t263 = qJD(2) * t245;
t259 = t243 * t266;
t233 = -t259 + t269;
t262 = qJD(4) * t233;
t261 = -r_i_i_C(3) - pkin(9) - pkin(2);
t260 = t243 * t269;
t258 = t242 * t265;
t257 = t242 * t264;
t256 = t242 * t263;
t255 = qJD(2) * t243 + qJD(1);
t254 = r_i_i_C(1) * t247 - r_i_i_C(2) * t244;
t253 = -t244 * r_i_i_C(1) - t247 * r_i_i_C(2);
t252 = qJ(3) - t253;
t250 = t254 * qJD(4) + qJD(3);
t232 = -qJD(1) * t260 - t246 * t263 + t255 * t266;
t230 = t251 * qJD(1) + t235 * qJD(2);
t229 = -qJD(1) * t259 - qJD(2) * t266 + t255 * t269;
t228 = t247 * t257 - t229 * t244 + (t235 * t247 - t244 * t272) * qJD(4);
t227 = -t244 * t257 - t229 * t247 + (-t235 * t244 - t247 * t272) * qJD(4);
t1 = [(-t247 * t262 - t274) * r_i_i_C(1) + (t244 * t262 - t273) * r_i_i_C(2) - t231 * qJ(3) - t233 * qJD(3) - pkin(1) * t264 + t261 * t232 + (t253 * t249 * qJD(4) + (-t254 - t275) * t265) * t242, t250 * (-t260 + t266) - t252 * t230 - t261 * t229, -t229, r_i_i_C(1) * t227 - t228 * r_i_i_C(2), 0, 0; t228 * r_i_i_C(1) + t227 * r_i_i_C(2) - t229 * qJ(3) + t235 * qJD(3) + t261 * t230 + (-pkin(1) * t246 + t275 * t270) * qJD(1), t261 * t231 + t252 * t232 + t250 * t251, t231 (-t244 * t258 + t273) * r_i_i_C(1) + (-t247 * t258 - t274) * r_i_i_C(2) + ((-t233 * t244 + t247 * t270) * r_i_i_C(1) + (-t233 * t247 - t244 * t270) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t250 * t245 + (t261 * t245 + t252 * t248) * qJD(2)) * t242, t256, t254 * t256 + ((-t243 * t247 + t244 * t271) * r_i_i_C(1) + (t243 * t244 + t247 * t271) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
