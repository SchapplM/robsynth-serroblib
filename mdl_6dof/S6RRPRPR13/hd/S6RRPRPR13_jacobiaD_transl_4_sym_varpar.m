% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR13
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR13_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR13_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR13_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR13_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:45:00
% EndTime: 2019-02-26 21:45:01
% DurationCPUTime: 0.18s
% Computational Cost: add. (179->58), mult. (534->98), div. (0->0), fcn. (502->8), ass. (0->43)
t274 = pkin(3) + pkin(8);
t242 = cos(pkin(6));
t245 = sin(qJ(1));
t247 = cos(qJ(2));
t266 = t245 * t247;
t244 = sin(qJ(2));
t248 = cos(qJ(1));
t267 = t244 * t248;
t234 = t242 * t266 + t267;
t250 = t242 * t267 + t266;
t230 = t234 * qJD(1) + t250 * qJD(2);
t243 = sin(qJ(4));
t273 = t230 * t243;
t246 = cos(qJ(4));
t272 = t230 * t246;
t241 = sin(pkin(6));
t271 = t241 * t245;
t270 = t241 * t247;
t269 = t241 * t248;
t268 = t244 * t245;
t265 = t247 * t248;
t264 = qJD(1) * t245;
t263 = qJD(1) * t248;
t262 = qJD(2) * t244;
t258 = t242 * t265;
t232 = -t258 + t268;
t261 = qJD(4) * t232;
t260 = -r_i_i_C(3) - pkin(9) - pkin(2);
t259 = t242 * t268;
t257 = t241 * t264;
t256 = t241 * t263;
t255 = t241 * t262;
t254 = qJD(2) * t242 + qJD(1);
t253 = r_i_i_C(1) * t246 - r_i_i_C(2) * t243;
t252 = -r_i_i_C(1) * t243 - r_i_i_C(2) * t246;
t251 = qJ(3) - t252;
t249 = t253 * qJD(4) + qJD(3);
t231 = -qJD(1) * t259 - t245 * t262 + t254 * t265;
t229 = t250 * qJD(1) + t234 * qJD(2);
t228 = -qJD(1) * t258 - qJD(2) * t265 + t254 * t268;
t227 = t246 * t256 - t228 * t243 + (t234 * t246 - t243 * t271) * qJD(4);
t226 = -t243 * t256 - t228 * t246 + (-t234 * t243 - t246 * t271) * qJD(4);
t1 = [(-t246 * t261 - t273) * r_i_i_C(1) + (t243 * t261 - t272) * r_i_i_C(2) - t230 * qJ(3) - t232 * qJD(3) - pkin(1) * t263 + t260 * t231 + (t252 * t248 * qJD(4) + (-t253 - t274) * t264) * t241, t249 * (-t259 + t265) - t251 * t229 - t260 * t228, -t228, r_i_i_C(1) * t226 - t227 * r_i_i_C(2), 0, 0; t227 * r_i_i_C(1) + t226 * r_i_i_C(2) - t228 * qJ(3) + t234 * qJD(3) + t260 * t229 + (-pkin(1) * t245 + t274 * t269) * qJD(1), t260 * t230 + t251 * t231 + t249 * t250, t230 (-t243 * t257 + t272) * r_i_i_C(1) + (-t246 * t257 - t273) * r_i_i_C(2) + ((-t232 * t243 + t246 * t269) * r_i_i_C(1) + (-t232 * t246 - t243 * t269) * r_i_i_C(2)) * qJD(4), 0, 0; 0 (t249 * t244 + (t260 * t244 + t251 * t247) * qJD(2)) * t241, t255, t253 * t255 + ((-t242 * t246 + t243 * t270) * r_i_i_C(1) + (t242 * t243 + t246 * t270) * r_i_i_C(2)) * qJD(4), 0, 0;];
JaD_transl  = t1;
