% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:01
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR2_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR2_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR2_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRRPR2_jacobiaD_transl_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:01:26
% EndTime: 2019-02-26 21:01:27
% DurationCPUTime: 0.21s
% Computational Cost: add. (297->54), mult. (374->82), div. (0->0), fcn. (292->10), ass. (0->46)
t218 = sin(qJ(3));
t219 = cos(qJ(4));
t249 = t219 * pkin(4);
t209 = pkin(3) + t249;
t214 = qJ(4) + pkin(11);
t210 = sin(t214);
t212 = cos(t214);
t233 = r_i_i_C(1) * t212 - r_i_i_C(2) * t210;
t229 = t209 + t233;
t220 = cos(qJ(3));
t248 = r_i_i_C(3) + qJ(5) + pkin(8);
t251 = t248 * t220;
t222 = -t229 * t218 + t251;
t256 = qJD(1) * t222;
t241 = t218 * qJD(5);
t242 = qJD(4) * t220;
t217 = sin(qJ(4));
t250 = pkin(4) * t217;
t255 = (-t209 * t218 + t251) * qJD(3) - t242 * t250 + t241;
t215 = qJ(1) + pkin(10);
t213 = cos(t215);
t234 = qJD(1) * t220 - qJD(4);
t231 = t234 * t213;
t235 = -qJD(1) + t242;
t252 = t235 * t212;
t246 = qJD(1) * t218;
t245 = qJD(3) * t218;
t244 = qJD(3) * t220;
t243 = qJD(4) * t218;
t240 = qJD(4) * t249;
t239 = pkin(7) + t250;
t238 = t248 * t218;
t232 = t235 * t210;
t211 = sin(t215);
t230 = t234 * t211;
t228 = r_i_i_C(1) * t210 + r_i_i_C(2) * t212 + t250;
t227 = -t209 * t220 - pkin(2) - t238;
t225 = t228 * t220;
t224 = t213 * t245 + t230;
t223 = t217 * t245 - t235 * t219;
t221 = qJD(5) * t220 + t228 * t243 + (-t229 * t220 - t238) * qJD(3);
t208 = -t212 * t231 + (t212 * t245 + t232) * t211;
t207 = t211 * t252 + (-t211 * t245 + t231) * t210;
t206 = t224 * t212 + t213 * t232;
t205 = t224 * t210 - t213 * t252;
t1 = [t213 * t240 + t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t255 * t211 + (-cos(qJ(1)) * pkin(1) - t239 * t211 + t227 * t213) * qJD(1), 0, -t211 * t256 + t221 * t213, t205 * r_i_i_C(1) + t206 * r_i_i_C(2) + (t223 * t213 + t217 * t230) * pkin(4), -t211 * t246 + t213 * t244, 0; t211 * t240 - t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t255 * t213 + (-sin(qJ(1)) * pkin(1) + t239 * t213 + t227 * t211) * qJD(1), 0, t221 * t211 + t213 * t256, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2) + (t223 * t211 - t217 * t231) * pkin(4), t211 * t244 + t213 * t246, 0; 0, 0, t222 * qJD(3) - qJD(4) * t225 + t241 (-t233 - t249) * t243 - qJD(3) * t225, t245, 0;];
JaD_transl  = t1;
