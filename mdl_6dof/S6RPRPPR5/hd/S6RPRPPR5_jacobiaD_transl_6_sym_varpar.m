% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:30
% EndTime: 2019-02-26 20:41:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (303->50), mult. (379->75), div. (0->0), fcn. (300->9), ass. (0->38)
t222 = pkin(9) + qJ(3);
t218 = sin(t222);
t220 = cos(t222);
t241 = pkin(5) * sin(pkin(10)) + qJ(4);
t245 = t220 * qJD(5);
t244 = pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
t254 = t244 * t218;
t261 = (-t220 * t241 + t254) * qJD(3) - (cos(pkin(10)) * pkin(5) + pkin(4) + pkin(7) + qJ(2)) * qJD(1) - t218 * qJD(4) - t245;
t221 = pkin(10) + qJ(6);
t217 = sin(t221);
t219 = cos(t221);
t234 = r_i_i_C(1) * t217 + r_i_i_C(2) * t219 + t241;
t259 = -t220 * t234 + t254;
t226 = sin(qJ(1));
t240 = qJD(6) * t218 + qJD(1);
t256 = t226 * t240;
t227 = cos(qJ(1));
t255 = t227 * t240;
t251 = qJD(1) * t226;
t250 = qJD(1) * t227;
t249 = qJD(3) * t218;
t248 = qJD(3) * t226;
t247 = qJD(3) * t227;
t246 = qJD(6) * t220;
t243 = t220 * t248;
t242 = t220 * t247;
t239 = -qJD(1) * t218 - qJD(6);
t237 = t244 * t220;
t233 = qJD(4) + (r_i_i_C(1) * t219 - r_i_i_C(2) * t217) * qJD(6);
t232 = t227 * t239 - t243;
t231 = t226 * t239 + t242;
t230 = qJD(2) + (-t241 * t218 - cos(pkin(9)) * pkin(2) - pkin(1) - t237) * qJD(1);
t228 = -qJD(5) * t218 + t233 * t220 + (-t218 * t234 - t237) * qJD(3);
t214 = t217 * t231 + t219 * t255;
t213 = -t217 * t255 + t219 * t231;
t212 = t217 * t232 - t219 * t256;
t211 = t217 * t256 + t219 * t232;
t1 = [t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t261 * t226 + t230 * t227, t250, t228 * t227 + t259 * t251, -t218 * t251 + t242, -t218 * t247 - t220 * t251, t213 * r_i_i_C(1) - t214 * r_i_i_C(2); t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t230 * t226 - t261 * t227, t251, t226 * t228 - t250 * t259, t218 * t250 + t243, -t218 * t248 + t220 * t250, -t211 * r_i_i_C(1) + t212 * r_i_i_C(2); 0, 0, -qJD(3) * t259 + t218 * t233 + t245, t249, qJD(3) * t220 (-t217 * t249 + t219 * t246) * r_i_i_C(2) + (t217 * t246 + t219 * t249) * r_i_i_C(1);];
JaD_transl  = t1;
