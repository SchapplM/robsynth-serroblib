% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6PRPRRR5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 19:56
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6PRPRRR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6PRPRRR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR5_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 19:56:11
% EndTime: 2019-02-26 19:56:11
% DurationCPUTime: 0.16s
% Computational Cost: add. (194->42), mult. (385->79), div. (0->0), fcn. (352->10), ass. (0->40)
t221 = qJ(4) + qJ(5);
t218 = sin(t221);
t220 = qJD(4) + qJD(5);
t253 = t218 * t220;
t219 = cos(t221);
t252 = t219 * t220;
t223 = sin(pkin(6));
t251 = t220 * t223;
t228 = cos(qJ(4));
t250 = t223 * t228;
t229 = cos(qJ(2));
t249 = t223 * t229;
t224 = cos(pkin(11));
t248 = t224 * t229;
t225 = cos(pkin(6));
t227 = sin(qJ(2));
t247 = t225 * t227;
t246 = t225 * t229;
t222 = sin(pkin(11));
t214 = t222 * t246 + t224 * t227;
t237 = qJD(2) * t248;
t242 = qJD(2) * t227;
t239 = t222 * t242;
t211 = -t225 * t239 + t237;
t235 = t222 * t251 - t211;
t245 = (-t214 * t253 - t235 * t219) * r_i_i_C(1) + (-t214 * t252 + t235 * t218) * r_i_i_C(2);
t212 = t222 * t227 - t224 * t246;
t234 = t222 * t229 + t224 * t247;
t209 = t234 * qJD(2);
t236 = t224 * t251 + t209;
t244 = (-t212 * t253 + t236 * t219) * r_i_i_C(1) + (-t212 * t252 - t236 * t218) * r_i_i_C(2);
t238 = t223 * t242;
t233 = -t220 * t225 + t238;
t240 = t220 * t249;
t243 = (t218 * t240 + t233 * t219) * r_i_i_C(1) + (-t233 * t218 + t219 * t240) * r_i_i_C(2);
t241 = -pkin(2) - r_i_i_C(3) - pkin(9) - pkin(8);
t226 = sin(qJ(4));
t232 = pkin(4) * t226 + r_i_i_C(1) * t218 + r_i_i_C(2) * t219 + qJ(3);
t231 = pkin(4) * qJD(4) * t228 + qJD(3) + (r_i_i_C(1) * t219 - r_i_i_C(2) * t218) * t220;
t1 = [0, t241 * t211 + t231 * (-t222 * t247 + t248) - t232 * t214 * qJD(2), t211 (t211 * t228 + (-t214 * t226 - t222 * t250) * qJD(4)) * pkin(4) + t245, t245, 0; 0, t241 * t209 + t231 * t234 - t232 * (-t225 * t237 + t239) t209 (t209 * t228 + (-t212 * t226 + t224 * t250) * qJD(4)) * pkin(4) + t244, t244, 0; 0 (t231 * t227 + (t241 * t227 + t232 * t229) * qJD(2)) * t223, t238 (t228 * t238 + (-t225 * t228 + t226 * t249) * qJD(4)) * pkin(4) + t243, t243, 0;];
JaD_transl  = t1;
