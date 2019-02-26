% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:47
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:57
% EndTime: 2019-02-26 20:46:57
% DurationCPUTime: 0.23s
% Computational Cost: add. (241->54), mult. (408->79), div. (0->0), fcn. (317->8), ass. (0->42)
t219 = cos(qJ(5));
t247 = t219 * pkin(5);
t209 = pkin(4) + t247;
t213 = qJ(3) + pkin(9);
t210 = sin(t213);
t211 = cos(t213);
t246 = r_i_i_C(3) + qJ(6) + pkin(8);
t228 = t246 * t211 - sin(qJ(3)) * pkin(3);
t263 = (-t209 * t210 - qJ(2) + t228) * qJD(1) - qJD(5) * t247 - qJD(4);
t216 = sin(qJ(5));
t249 = r_i_i_C(2) * t216;
t230 = r_i_i_C(1) * t219 + t209 - t249;
t262 = -(-t230 * t210 + t228) * qJD(3) - qJD(6) * t210;
t218 = sin(qJ(1));
t233 = qJD(1) * t210 + qJD(5);
t221 = cos(qJ(1));
t240 = qJD(3) * t221;
t254 = -t211 * t240 + t233 * t218;
t261 = t254 * t216;
t226 = t246 * t210 + cos(qJ(3)) * pkin(3);
t258 = t230 * t211 + t226;
t253 = pkin(5) + r_i_i_C(1);
t229 = r_i_i_C(2) * t219 + t253 * t216;
t250 = pkin(5) * t216;
t244 = t219 * t221;
t243 = qJD(1) * t218;
t212 = qJD(1) * t221;
t242 = qJD(3) * t211;
t241 = qJD(3) * t218;
t239 = qJD(5) * t210;
t238 = qJD(5) * t211;
t236 = t211 * qJD(6);
t234 = qJD(1) + t239;
t231 = t234 * t221;
t225 = qJD(5) * t229;
t223 = qJD(1) * t258;
t207 = -t234 * t219 * t218 + (-t211 * t241 - t233 * t221) * t216;
t222 = -t239 * t250 - t236 + qJD(2) + (-pkin(1) - qJ(4) - pkin(7) - t250) * qJD(1) + (t209 * t211 + t226) * qJD(3);
t208 = t233 * t244 + (-t234 * t216 + t219 * t242) * t218;
t206 = -t216 * t231 - t219 * t254;
t205 = -t219 * t231 + t261;
t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t263 * t218 + t222 * t221, t212, t221 * t223 + (-t229 * t238 - t262) * t218, -t243, -t208 * r_i_i_C(2) + t253 * t207, t210 * t241 - t211 * t212; t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t222 * t218 - t263 * t221, t243, t218 * t223 + (t211 * t225 + t262) * t221, t212, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t234 * t244 - t261) * pkin(5), -t210 * t240 - t211 * t243; 0, 0, -t258 * qJD(3) + t210 * t225 + t236, 0 (-t253 * t219 + t249) * t238 + t229 * t210 * qJD(3), t242;];
JaD_transl  = t1;
