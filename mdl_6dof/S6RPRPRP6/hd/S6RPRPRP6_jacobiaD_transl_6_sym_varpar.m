% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:28
% EndTime: 2019-02-26 20:46:28
% DurationCPUTime: 0.26s
% Computational Cost: add. (275->53), mult. (446->78), div. (0->0), fcn. (348->7), ass. (0->42)
t214 = sin(qJ(1));
t215 = cos(qJ(5));
t213 = sin(qJ(5));
t210 = pkin(9) + qJ(3);
t208 = sin(t210);
t232 = qJD(5) * t208 + qJD(1);
t226 = t232 * t213;
t209 = cos(t210);
t241 = qJD(3) * t209;
t257 = (-t215 * t241 + t226) * t214;
t249 = pkin(5) + r_i_i_C(1);
t223 = r_i_i_C(2) * t215 + t249 * t213;
t251 = qJ(4) + t223;
t237 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t252 = t237 * t208;
t219 = t251 * t209 - t252;
t233 = pkin(5) * t213 + qJ(4);
t238 = t209 * qJD(6);
t247 = pkin(5) * qJD(5);
t256 = (-t233 * t209 + t252) * qJD(3) - (t215 * pkin(5) + pkin(4) + pkin(7) + qJ(2)) * qJD(1) - (t215 * t247 + qJD(4)) * t208 - t238;
t216 = cos(qJ(1));
t246 = t215 * t216;
t244 = qJD(1) * t214;
t243 = qJD(1) * t216;
t242 = qJD(3) * t208;
t240 = qJD(3) * t214;
t239 = qJD(3) * t216;
t235 = t209 * t240;
t234 = t209 * t239;
t231 = qJD(1) * t208 + qJD(5);
t229 = t237 * t209;
t225 = t231 * t216;
t224 = -r_i_i_C(2) * t213 + t249 * t215;
t221 = t224 * qJD(5) + qJD(4);
t220 = -t231 * t214 + t234;
t218 = -t213 * t247 + qJD(2) + (-t233 * t208 - cos(pkin(9)) * pkin(2) - pkin(1) - t229) * qJD(1);
t204 = t220 * t215 - t216 * t226;
t217 = -qJD(6) * t208 + t221 * t209 + (-t208 * t251 - t229) * qJD(3);
t205 = t220 * t213 + t232 * t246;
t203 = -t232 * t215 * t214 + (-t225 - t235) * t213;
t202 = -t215 * t225 + t257;
t1 = [t203 * r_i_i_C(1) + t202 * r_i_i_C(2) + t256 * t214 + t218 * t216, t243, t217 * t216 - t219 * t244, -t208 * t244 + t234, -t205 * r_i_i_C(2) + t249 * t204, -t208 * t239 - t209 * t244; t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t218 * t214 - t256 * t216, t244, t217 * t214 + t219 * t243, t208 * t243 + t235, -t202 * r_i_i_C(1) + t203 * r_i_i_C(2) + (t231 * t246 - t257) * pkin(5), -t208 * t240 + t209 * t243; 0, 0, t219 * qJD(3) + t221 * t208 + t238, t242, t223 * t209 * qJD(5) + t224 * t242, t241;];
JaD_transl  = t1;
