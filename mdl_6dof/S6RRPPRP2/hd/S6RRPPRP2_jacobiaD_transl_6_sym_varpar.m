% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
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

function JaD_transl = S6RRPPRP2_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP2_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP2_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPRP2_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:25:22
% EndTime: 2019-02-26 21:25:23
% DurationCPUTime: 0.28s
% Computational Cost: add. (282->54), mult. (466->79), div. (0->0), fcn. (361->8), ass. (0->42)
t213 = qJ(2) + pkin(9);
t211 = sin(t213);
t212 = cos(t213);
t219 = cos(qJ(5));
t242 = pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t228 = t242 * t211 + sin(qJ(2)) * pkin(2);
t216 = sin(qJ(5));
t238 = pkin(5) * t216 + qJ(4);
t243 = t212 * qJD(6);
t252 = pkin(5) * qJD(5);
t264 = (-t238 * t212 + t228) * qJD(2) - (t219 * pkin(5) + pkin(4) + pkin(7) + qJ(3)) * qJD(1) - (t219 * t252 + qJD(4)) * t211 - t243;
t218 = sin(qJ(1));
t237 = qJD(5) * t211 + qJD(1);
t233 = t237 * t216;
t246 = qJD(2) * t212;
t263 = (-t219 * t246 + t233) * t218;
t256 = pkin(5) + r_i_i_C(1);
t230 = r_i_i_C(2) * t219 + t256 * t216;
t259 = qJ(4) + t230;
t224 = t259 * t212 - t228;
t258 = -t242 * t212 - cos(qJ(2)) * pkin(2);
t221 = cos(qJ(1));
t251 = t219 * t221;
t249 = qJD(1) * t218;
t248 = qJD(1) * t221;
t247 = qJD(2) * t211;
t245 = qJD(2) * t218;
t244 = qJD(2) * t221;
t240 = t212 * t245;
t239 = t212 * t244;
t236 = qJD(1) * t211 + qJD(5);
t232 = t236 * t221;
t231 = -r_i_i_C(2) * t216 + t256 * t219;
t226 = t231 * qJD(5) + qJD(4);
t225 = -t236 * t218 + t239;
t223 = -t216 * t252 + qJD(3) + (-t238 * t211 - pkin(1) + t258) * qJD(1);
t207 = t225 * t219 - t221 * t233;
t222 = -qJD(6) * t211 + t226 * t212 + (-t211 * t259 + t258) * qJD(2);
t208 = t225 * t216 + t237 * t251;
t206 = -t237 * t219 * t218 + (-t232 - t240) * t216;
t205 = -t219 * t232 + t263;
t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t264 * t218 + t223 * t221, t222 * t221 - t224 * t249, t248, -t211 * t249 + t239, -t208 * r_i_i_C(2) + t256 * t207, -t211 * t244 - t212 * t249; t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t223 * t218 - t264 * t221, t222 * t218 + t224 * t248, t249, t211 * t248 + t240, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t236 * t251 - t263) * pkin(5), -t211 * t245 + t212 * t248; 0, t224 * qJD(2) + t226 * t211 + t243, 0, t247, t230 * t212 * qJD(5) + t231 * t247, t246;];
JaD_transl  = t1;
