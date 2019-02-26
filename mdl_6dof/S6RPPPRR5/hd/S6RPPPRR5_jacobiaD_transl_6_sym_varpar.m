% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:24:54
% EndTime: 2019-02-26 20:24:54
% DurationCPUTime: 0.26s
% Computational Cost: add. (176->54), mult. (490->85), div. (0->0), fcn. (467->8), ass. (0->40)
t224 = cos(qJ(5));
t221 = sin(qJ(6));
t223 = cos(qJ(6));
t232 = t223 * r_i_i_C(1) - t221 * r_i_i_C(2) + pkin(5);
t222 = sin(qJ(5));
t248 = pkin(8) + r_i_i_C(3);
t239 = t248 * t222;
t252 = t232 * t224 + t239;
t253 = t252 * qJD(5);
t251 = -pkin(1) - qJ(3);
t250 = pkin(3) + qJ(2);
t247 = sin(qJ(1));
t246 = sin(pkin(9));
t225 = cos(qJ(1));
t219 = qJD(1) * t225;
t245 = qJD(5) * t222;
t244 = qJD(5) * t224;
t220 = cos(pkin(9));
t216 = t247 * t220 + t225 * t246;
t243 = qJD(6) * t216;
t242 = qJD(6) * t222;
t241 = qJD(6) * t223;
t240 = qJD(6) * t224;
t238 = t248 * t224;
t237 = qJD(1) * t247;
t236 = t247 * t246;
t235 = t221 * r_i_i_C(1) + t223 * r_i_i_C(2);
t213 = t216 * qJD(1);
t234 = -t216 * t240 + t213;
t214 = -qJD(1) * t236 + t220 * t219;
t215 = t225 * t220 - t236;
t233 = t215 * t240 - t214;
t231 = qJD(6) * t235;
t230 = -t213 * t224 - t215 * t245 + t243;
t229 = qJD(6) * t215 - t214 * t224 + t216 * t245;
t227 = -t232 * t222 + t238;
t226 = t227 * qJD(5) - t224 * t231;
t212 = t234 * t221 - t229 * t223;
t211 = t229 * t221 + t234 * t223;
t1 = [(t214 * t221 + t216 * t241) * r_i_i_C(1) + (t214 * t223 - t221 * t243) * r_i_i_C(2) + t214 * pkin(7) - t247 * qJD(3) + t225 * qJD(2) + (-pkin(4) - t252) * t213 + (t251 * t225 - t250 * t247) * qJD(1) + t226 * t215, t219, -t237, 0, t227 * t214 + (t235 * t242 - t253) * t216, t211 * r_i_i_C(1) - t212 * r_i_i_C(2); t247 * qJD(2) + t213 * pkin(7) + t212 * r_i_i_C(1) + t211 * r_i_i_C(2) + t225 * qJD(3) + (-pkin(5) * t222 + t238) * t216 * qJD(5) + (pkin(5) * t224 + pkin(4) + t239) * t214 + (t250 * t225 + t251 * t247) * qJD(1), t237, t219, 0, t227 * t213 + (-t222 * t231 + t253) * t215 (t233 * r_i_i_C(1) + t230 * r_i_i_C(2)) * t223 + (t230 * r_i_i_C(1) - t233 * r_i_i_C(2)) * t221; 0, 0, 0, 0, t226 (t221 * t242 - t223 * t244) * r_i_i_C(2) + (-t221 * t244 - t222 * t241) * r_i_i_C(1);];
JaD_transl  = t1;
