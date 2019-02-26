% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiaD_transl_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:12
% EndTime: 2019-02-26 20:39:12
% DurationCPUTime: 0.20s
% Computational Cost: add. (366->52), mult. (333->79), div. (0->0), fcn. (263->12), ass. (0->38)
t221 = cos(pkin(11)) * pkin(5) + pkin(4);
t230 = qJ(3) + pkin(10);
t224 = sin(t230);
t227 = cos(t230);
t260 = r_i_i_C(3) + pkin(8) + qJ(5);
t242 = t260 * t227 - sin(qJ(3)) * pkin(3);
t250 = t224 * qJD(5);
t266 = t250 + (-t221 * t224 + t242) * qJD(3);
t229 = pkin(11) + qJ(6);
t223 = sin(t229);
t226 = cos(t229);
t243 = r_i_i_C(1) * t226 - r_i_i_C(2) * t223 + t221;
t238 = -t243 * t224 + t242;
t263 = -t260 * t224 - cos(qJ(3)) * pkin(3);
t231 = qJ(1) + pkin(9);
t228 = cos(t231);
t258 = t226 * t228;
t225 = sin(t231);
t257 = qJD(1) * t225;
t256 = qJD(1) * t228;
t255 = qJD(3) * t224;
t254 = qJD(3) * t227;
t253 = qJD(3) * t228;
t252 = qJD(6) * t224;
t251 = qJD(6) * t227;
t248 = pkin(5) * sin(pkin(11)) + qJ(4) + pkin(7);
t247 = -qJD(1) + t251;
t246 = qJD(1) * t227 - qJD(6);
t245 = r_i_i_C(1) * t223 + r_i_i_C(2) * t226;
t244 = t247 * t223;
t240 = -t221 * t227 - pkin(2) + t263;
t239 = t224 * t253 + t246 * t225;
t237 = qJD(5) * t227 + t245 * t252 + (-t243 * t227 + t263) * qJD(3);
t220 = -t246 * t258 + (t226 * t255 + t244) * t225;
t219 = t247 * t226 * t225 + (-t225 * t255 + t246 * t228) * t223;
t218 = t239 * t226 + t228 * t244;
t217 = t239 * t223 - t247 * t258;
t1 = [t220 * r_i_i_C(1) + t219 * r_i_i_C(2) + t228 * qJD(4) - t266 * t225 + (-cos(qJ(1)) * pkin(1) - t248 * t225 + t240 * t228) * qJD(1), 0, t237 * t228 - t238 * t257, t256, -t224 * t257 + t227 * t253, t217 * r_i_i_C(1) + t218 * r_i_i_C(2); -t218 * r_i_i_C(1) + t217 * r_i_i_C(2) + t225 * qJD(4) + t266 * t228 + (-sin(qJ(1)) * pkin(1) + t248 * t228 + t240 * t225) * qJD(1), 0, t237 * t225 + t238 * t256, t257, t224 * t256 + t225 * t254, -t219 * r_i_i_C(1) + t220 * r_i_i_C(2); 0, 0, t238 * qJD(3) - t245 * t251 + t250, 0, t255 (t223 * t252 - t226 * t254) * r_i_i_C(2) + (-t223 * t254 - t226 * t252) * r_i_i_C(1);];
JaD_transl  = t1;
