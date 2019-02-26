% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:23:11
% EndTime: 2019-02-26 20:23:11
% DurationCPUTime: 0.15s
% Computational Cost: add. (259->47), mult. (292->74), div. (0->0), fcn. (227->9), ass. (0->34)
t217 = pkin(10) + qJ(5);
t213 = sin(t217);
t215 = cos(t217);
t239 = pkin(8) + r_i_i_C(3);
t232 = t239 * t215;
t244 = t232 - pkin(4) * sin(pkin(10)) - pkin(5) * t213 - qJ(3);
t221 = sin(qJ(6));
t222 = cos(qJ(6));
t228 = r_i_i_C(1) * t222 - r_i_i_C(2) * t221 + pkin(5);
t233 = t239 * t213;
t243 = t228 * t215 + t233;
t230 = qJD(1) * t213 + qJD(6);
t242 = t221 * t230;
t241 = t222 * t230;
t238 = -pkin(2) - pkin(7) - qJ(4);
t218 = qJ(1) + pkin(9);
t214 = sin(t218);
t237 = qJD(1) * t214;
t216 = cos(t218);
t212 = qJD(1) * t216;
t236 = qJD(5) * t221;
t235 = qJD(5) * t222;
t234 = qJD(6) * t215;
t231 = -qJD(6) * t213 - qJD(1);
t229 = r_i_i_C(1) * t221 + r_i_i_C(2) * t222;
t226 = t229 * qJD(6);
t225 = qJD(3) + (pkin(5) * t215 + t233) * qJD(5);
t224 = -t215 * t236 + t231 * t222;
t223 = t215 * t235 + t231 * t221;
t211 = t223 * t214 + t216 * t241;
t210 = t224 * t214 - t216 * t242;
t209 = -t214 * t241 + t223 * t216;
t208 = t214 * t242 + t224 * t216;
t1 = [t209 * r_i_i_C(1) + t208 * r_i_i_C(2) - t214 * qJD(4) + t225 * t216 + (-cos(qJ(1)) * pkin(1) + t238 * t216 + t244 * t214) * qJD(1), 0, t212, -t237, t243 * t212 + (-t229 * t234 + (-t228 * t213 + t232) * qJD(5)) * t214, t210 * r_i_i_C(1) - t211 * r_i_i_C(2); t211 * r_i_i_C(1) + t210 * r_i_i_C(2) + t216 * qJD(4) + t225 * t214 + (-sin(qJ(1)) * pkin(1) + t238 * t214 - t244 * t216) * qJD(1), 0, t237, t212 (t228 * t216 * qJD(5) + t239 * t237) * t213 + (t228 * t237 + (-t239 * qJD(5) + t226) * t216) * t215, -t208 * r_i_i_C(1) + t209 * r_i_i_C(2); 0, 0, 0, 0, -t243 * qJD(5) + t213 * t226 (t213 * t235 + t221 * t234) * r_i_i_C(2) + (t213 * t236 - t222 * t234) * r_i_i_C(1);];
JaD_transl  = t1;
