% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:25
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:22
% EndTime: 2019-02-26 21:25:22
% DurationCPUTime: 0.18s
% Computational Cost: add. (204->44), mult. (344->68), div. (0->0), fcn. (267->8), ass. (0->35)
t208 = qJ(2) + pkin(9);
t206 = sin(t208);
t207 = cos(t208);
t230 = pkin(3) + pkin(8) + r_i_i_C(3);
t222 = t230 * t206 + sin(qJ(2)) * pkin(2);
t249 = (-qJ(4) * t207 + t222) * qJD(2) - (pkin(4) + qJ(3) + pkin(7)) * qJD(1) - t206 * qJD(4);
t210 = sin(qJ(5));
t213 = cos(qJ(5));
t223 = r_i_i_C(1) * t210 + r_i_i_C(2) * t213 + qJ(4);
t247 = -t223 * t207 + t222;
t226 = qJD(5) * t206 + qJD(1);
t246 = t213 * t226;
t244 = t226 * t210;
t243 = -t230 * t207 - cos(qJ(2)) * pkin(2);
t212 = sin(qJ(1));
t237 = qJD(1) * t212;
t215 = cos(qJ(1));
t236 = qJD(1) * t215;
t235 = qJD(2) * t206;
t234 = qJD(2) * t207;
t233 = qJD(2) * t213;
t232 = qJD(5) * t207;
t229 = t212 * t234;
t228 = t215 * t234;
t225 = -qJD(1) * t206 - qJD(5);
t224 = t225 * t215;
t220 = qJD(4) + (r_i_i_C(1) * t213 - r_i_i_C(2) * t210) * qJD(5);
t219 = t225 * t212 + t228;
t218 = qJD(3) + (-qJ(4) * t206 - pkin(1) + t243) * qJD(1);
t216 = t220 * t207 + (-t223 * t206 + t243) * qJD(2);
t204 = t219 * t210 + t215 * t246;
t203 = t219 * t213 - t215 * t244;
t202 = -t212 * t246 + (t224 - t229) * t210;
t201 = t213 * t224 + (-t207 * t233 + t244) * t212;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t249 * t212 + t218 * t215, t216 * t215 + t247 * t237, t236, -t206 * t237 + t228, t203 * r_i_i_C(1) - t204 * r_i_i_C(2), 0; t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t218 * t212 - t249 * t215, t216 * t212 - t236 * t247, t237, t206 * t236 + t229, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0; 0, -qJD(2) * t247 + t220 * t206, 0, t235 (-t210 * t235 + t213 * t232) * r_i_i_C(2) + (t206 * t233 + t210 * t232) * r_i_i_C(1), 0;];
JaD_transl  = t1;
