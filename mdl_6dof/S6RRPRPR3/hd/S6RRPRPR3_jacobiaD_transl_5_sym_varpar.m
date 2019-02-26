% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRPR3_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:38:52
% EndTime: 2019-02-26 21:38:52
% DurationCPUTime: 0.26s
% Computational Cost: add. (296->51), mult. (396->76), div. (0->0), fcn. (309->10), ass. (0->41)
t230 = cos(qJ(4));
t260 = pkin(4) * t230;
t217 = pkin(3) + t260;
t224 = qJ(2) + pkin(10);
t220 = sin(t224);
t222 = cos(t224);
t258 = r_i_i_C(3) + qJ(5) + pkin(8);
t241 = t258 * t222 - sin(qJ(2)) * pkin(2);
t248 = qJD(4) * t222 - qJD(1);
t251 = t220 * qJD(5);
t227 = sin(qJ(4));
t261 = pkin(4) * t227;
t268 = (-t217 * t220 + t241) * qJD(2) - t248 * t261 - qJD(1) * (-qJ(3) - pkin(7)) + t251;
t223 = qJ(4) + pkin(11);
t219 = sin(t223);
t221 = cos(t223);
t246 = r_i_i_C(1) * t221 - r_i_i_C(2) * t219;
t242 = t217 + t246;
t235 = -t242 * t220 + t241;
t229 = sin(qJ(1));
t244 = t248 * t229;
t232 = cos(qJ(1));
t245 = t248 * t232;
t265 = -t258 * t220 - cos(qJ(2)) * pkin(2);
t247 = qJD(1) * t222 - qJD(4);
t254 = qJD(2) * t229;
t264 = -t220 * t254 + t247 * t232;
t256 = qJD(1) * t229;
t255 = qJD(1) * t232;
t253 = qJD(2) * t232;
t252 = qJD(4) * t220;
t239 = r_i_i_C(1) * t219 + r_i_i_C(2) * t221 + t261;
t238 = t239 * t222;
t236 = t220 * t253 + t247 * t229;
t234 = qJD(4) * t260 + qJD(3) + (-t217 * t222 - pkin(1) + t265) * qJD(1);
t233 = qJD(5) * t222 + t239 * t252 + (-t242 * t222 + t265) * qJD(2);
t216 = t219 * t244 - t221 * t264;
t215 = t264 * t219 + t221 * t244;
t214 = t219 * t245 + t236 * t221;
t213 = t236 * t219 - t221 * t245;
t1 = [t216 * r_i_i_C(1) + t215 * r_i_i_C(2) - t268 * t229 + t234 * t232, t233 * t232 - t235 * t256, t255, t213 * r_i_i_C(1) + t214 * r_i_i_C(2) + (t236 * t227 - t230 * t245) * pkin(4), -t220 * t256 + t222 * t253, 0; -t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t234 * t229 + t268 * t232, t233 * t229 + t235 * t255, t256, -t215 * r_i_i_C(1) + t216 * r_i_i_C(2) + (-t227 * t264 - t230 * t244) * pkin(4), t220 * t255 + t222 * t254, 0; 0, t235 * qJD(2) - qJD(4) * t238 + t251, 0 (-t246 - t260) * t252 - qJD(2) * t238, qJD(2) * t220, 0;];
JaD_transl  = t1;
