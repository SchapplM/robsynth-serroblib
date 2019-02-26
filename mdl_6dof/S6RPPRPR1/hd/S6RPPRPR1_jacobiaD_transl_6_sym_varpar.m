% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:25:30
% EndTime: 2019-02-26 20:25:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (353->51), mult. (313->78), div. (0->0), fcn. (250->11), ass. (0->38)
t218 = cos(pkin(11)) * pkin(5) + pkin(4);
t227 = pkin(10) + qJ(4);
t221 = sin(t227);
t245 = t221 * qJD(5);
t224 = cos(t227);
t255 = r_i_i_C(3) + pkin(8) + qJ(5);
t256 = t255 * t224;
t259 = (-t218 * t221 + t256) * qJD(4) + t245;
t226 = pkin(11) + qJ(6);
t220 = sin(t226);
t223 = cos(t226);
t236 = r_i_i_C(1) * t223 - r_i_i_C(2) * t220 + t218;
t233 = -t236 * t221 + t256;
t228 = qJ(1) + pkin(9);
t225 = cos(t228);
t253 = t223 * t225;
t222 = sin(t228);
t252 = qJD(1) * t222;
t251 = qJD(1) * t225;
t250 = qJD(4) * t221;
t249 = qJD(4) * t224;
t248 = qJD(4) * t225;
t247 = qJD(6) * t221;
t246 = qJD(6) * t224;
t244 = t255 * t221;
t241 = pkin(5) * sin(pkin(11)) + pkin(7) + qJ(3);
t240 = -qJD(1) + t246;
t239 = qJD(1) * t224 - qJD(6);
t238 = r_i_i_C(1) * t220 + r_i_i_C(2) * t223;
t237 = t240 * t220;
t235 = -t218 * t224 - cos(pkin(10)) * pkin(3) - pkin(2) - t244;
t234 = t221 * t248 + t239 * t222;
t232 = qJD(5) * t224 + t238 * t247 + (-t236 * t224 - t244) * qJD(4);
t217 = -t239 * t253 + (t223 * t250 + t237) * t222;
t216 = t240 * t223 * t222 + (-t222 * t250 + t239 * t225) * t220;
t215 = t234 * t223 + t225 * t237;
t214 = t234 * t220 - t240 * t253;
t1 = [t217 * r_i_i_C(1) + t216 * r_i_i_C(2) + t225 * qJD(3) - t259 * t222 + (-cos(qJ(1)) * pkin(1) - t241 * t222 + t235 * t225) * qJD(1), 0, t251, t232 * t225 - t233 * t252, -t221 * t252 + t224 * t248, t214 * r_i_i_C(1) + t215 * r_i_i_C(2); -t215 * r_i_i_C(1) + t214 * r_i_i_C(2) + t222 * qJD(3) + t259 * t225 + (-sin(qJ(1)) * pkin(1) + t241 * t225 + t235 * t222) * qJD(1), 0, t252, t232 * t222 + t233 * t251, t221 * t251 + t222 * t249, -t216 * r_i_i_C(1) + t217 * r_i_i_C(2); 0, 0, 0, t233 * qJD(4) - t238 * t246 + t245, t250 (t220 * t247 - t223 * t249) * r_i_i_C(2) + (-t220 * t249 - t223 * t247) * r_i_i_C(1);];
JaD_transl  = t1;
