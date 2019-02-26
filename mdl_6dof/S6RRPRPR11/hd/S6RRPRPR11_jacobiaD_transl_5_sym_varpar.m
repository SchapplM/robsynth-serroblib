% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPRPR11
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:43
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR11_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR11_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR11_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR11_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:43:37
% EndTime: 2019-02-26 21:43:37
% DurationCPUTime: 0.26s
% Computational Cost: add. (207->53), mult. (440->80), div. (0->0), fcn. (342->8), ass. (0->44)
t209 = sin(qJ(2));
t211 = cos(qJ(4));
t212 = cos(qJ(2));
t208 = sin(qJ(4));
t243 = pkin(4) * t208;
t230 = qJ(3) + t243;
t234 = t212 * qJD(5);
t240 = pkin(4) * qJD(4);
t242 = t211 * pkin(4);
t233 = pkin(2) + r_i_i_C(3) + qJ(5) + pkin(8);
t245 = t233 * t209;
t252 = (-t230 * t212 + t245) * qJD(2) - (pkin(7) + pkin(3) + t242) * qJD(1) - (t211 * t240 + qJD(3)) * t209 - t234;
t206 = qJ(4) + pkin(10);
t204 = sin(t206);
t205 = cos(t206);
t222 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205 + t243;
t220 = qJ(3) + t222;
t250 = -t220 * t212 + t245;
t210 = sin(qJ(1));
t229 = qJD(4) * t209 + qJD(1);
t249 = t210 * t229;
t213 = cos(qJ(1));
t228 = qJD(1) * t209 + qJD(4);
t247 = t228 * t213;
t223 = t229 * t213;
t239 = qJD(1) * t210;
t238 = qJD(1) * t213;
t237 = qJD(2) * t209;
t236 = qJD(2) * t212;
t235 = qJD(2) * t213;
t232 = t210 * t236;
t231 = t212 * t235;
t226 = t233 * t212;
t221 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t242;
t219 = -t232 - t247;
t218 = -t228 * t210 + t231;
t217 = t221 * qJD(4) + qJD(3);
t215 = -t208 * t240 + (-t230 * t209 - pkin(1) - t226) * qJD(1);
t214 = -qJD(5) * t209 + t217 * t212 + (-t220 * t209 - t226) * qJD(2);
t202 = t218 * t204 + t205 * t223;
t201 = -t204 * t223 + t218 * t205;
t200 = t219 * t204 - t205 * t249;
t199 = t204 * t249 + t219 * t205;
t1 = [t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t252 * t210 + t215 * t213, t214 * t213 + t250 * t239, -t209 * t239 + t231, t201 * r_i_i_C(1) - t202 * r_i_i_C(2) + (-t208 * t223 + t218 * t211) * pkin(4), -t209 * t235 - t212 * t239, 0; t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t215 * t210 - t252 * t213, t214 * t210 - t238 * t250, t209 * t238 + t232, -t199 * r_i_i_C(1) + t200 * r_i_i_C(2) + (t211 * t247 + (-t229 * t208 + t211 * t236) * t210) * pkin(4), -t210 * t237 + t212 * t238, 0; 0, -qJD(2) * t250 + t217 * t209 + t234, t237, t222 * t212 * qJD(4) + t221 * t237, t236, 0;];
JaD_transl  = t1;
