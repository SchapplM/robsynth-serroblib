% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:01
% EndTime: 2019-02-26 20:28:02
% DurationCPUTime: 0.18s
% Computational Cost: add. (179->48), mult. (319->73), div. (0->0), fcn. (254->8), ass. (0->37)
t203 = cos(pkin(9)) * pkin(5) + pkin(4);
t212 = cos(qJ(4));
t229 = t212 * qJD(5);
t210 = sin(qJ(4));
t237 = r_i_i_C(3) + pkin(8) + qJ(5);
t239 = t237 * t210;
t244 = (t203 * t212 + t239) * qJD(4) + qJD(3) - t229;
t207 = pkin(9) + qJ(6);
t204 = sin(t207);
t205 = cos(t207);
t218 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
t242 = t218 * t212 + t239;
t213 = cos(qJ(1));
t231 = qJD(6) * t210;
t224 = qJD(1) + t231;
t241 = t213 * t224;
t211 = sin(qJ(1));
t223 = qJD(1) * t210 + qJD(6);
t233 = qJD(4) * t212;
t238 = t211 * t233 + t223 * t213;
t235 = qJD(1) * t211;
t206 = qJD(1) * t213;
t234 = qJD(4) * t210;
t232 = qJD(4) * t213;
t230 = qJD(6) * t212;
t225 = t237 * t212;
t221 = pkin(5) * sin(pkin(9)) + pkin(7) - qJ(2);
t220 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
t219 = t224 * t211;
t217 = -t203 * t210 - pkin(1) - qJ(3) + t225;
t216 = t223 * t211 - t212 * t232;
t214 = qJD(5) * t210 - t220 * t230 + (-t218 * t210 + t225) * qJD(4);
t202 = t204 * t219 - t238 * t205;
t201 = t238 * t204 + t205 * t219;
t200 = t204 * t241 + t216 * t205;
t199 = t216 * t204 - t205 * t241;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t213 * qJD(2) - t244 * t211 + (t221 * t211 + t217 * t213) * qJD(1), t206, -t235, t214 * t213 - t235 * t242, t210 * t232 + t212 * t235, t199 * r_i_i_C(1) + t200 * r_i_i_C(2); -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t211 * qJD(2) + t244 * t213 + (t217 * t211 - t221 * t213) * qJD(1), t235, t206, t242 * t206 + t214 * t211, -t212 * t206 + t211 * t234, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2); 0, 0, 0, -qJD(4) * t242 + t220 * t231 + t229, t233 (t204 * t230 + t205 * t234) * r_i_i_C(2) + (t204 * t234 - t205 * t230) * r_i_i_C(1);];
JaD_transl  = t1;
