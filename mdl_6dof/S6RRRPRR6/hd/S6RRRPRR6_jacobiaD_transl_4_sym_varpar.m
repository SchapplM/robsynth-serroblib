% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:18
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRPRR6_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRR6_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRPRR6_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPRR6_jacobiaD_transl_4_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:18:48
% EndTime: 2019-02-26 22:18:49
% DurationCPUTime: 0.19s
% Computational Cost: add. (189->46), mult. (370->74), div. (0->0), fcn. (290->8), ass. (0->40)
t213 = cos(qJ(3));
t242 = t213 * pkin(3);
t205 = pkin(2) + t242;
t211 = sin(qJ(2));
t214 = cos(qJ(2));
t229 = qJD(3) * t214 - qJD(1);
t234 = t211 * qJD(4);
t210 = sin(qJ(3));
t243 = pkin(3) * t210;
t241 = r_i_i_C(3) + qJ(4) + pkin(8);
t246 = t241 * t214;
t249 = pkin(7) * qJD(1) + (-t205 * t211 + t246) * qJD(2) - t229 * t243 + t234;
t208 = qJ(3) + pkin(11);
t206 = sin(t208);
t207 = cos(t208);
t227 = r_i_i_C(1) * t207 - r_i_i_C(2) * t206;
t223 = t205 + t227;
t218 = -t223 * t211 + t246;
t212 = sin(qJ(1));
t224 = t229 * t212;
t215 = cos(qJ(1));
t225 = t229 * t215;
t228 = qJD(1) * t214 - qJD(3);
t237 = qJD(2) * t212;
t245 = -t211 * t237 + t228 * t215;
t239 = qJD(1) * t212;
t238 = qJD(1) * t215;
t236 = qJD(2) * t215;
t235 = qJD(3) * t211;
t232 = t241 * t211;
t222 = r_i_i_C(1) * t206 + r_i_i_C(2) * t207 + t243;
t221 = t222 * t214;
t219 = t211 * t236 + t228 * t212;
t217 = qJD(3) * t242 + (-t205 * t214 - pkin(1) - t232) * qJD(1);
t216 = qJD(4) * t214 + t222 * t235 + (-t223 * t214 - t232) * qJD(2);
t204 = t206 * t224 - t207 * t245;
t203 = t245 * t206 + t207 * t224;
t202 = t206 * t225 + t219 * t207;
t201 = t219 * t206 - t207 * t225;
t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t249 * t212 + t217 * t215, t216 * t215 - t218 * t239, t201 * r_i_i_C(1) + t202 * r_i_i_C(2) + (t219 * t210 - t213 * t225) * pkin(3), -t211 * t239 + t214 * t236, 0, 0; -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t217 * t212 + t249 * t215, t216 * t212 + t218 * t238, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2) + (-t210 * t245 - t213 * t224) * pkin(3), t211 * t238 + t214 * t237, 0, 0; 0, t218 * qJD(2) - qJD(3) * t221 + t234 (-t227 - t242) * t235 - qJD(2) * t221, qJD(2) * t211, 0, 0;];
JaD_transl  = t1;
