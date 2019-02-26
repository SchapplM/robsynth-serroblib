% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:33
% EndTime: 2019-02-26 20:32:33
% DurationCPUTime: 0.20s
% Computational Cost: add. (138->48), mult. (386->70), div. (0->0), fcn. (302->6), ass. (0->39)
t241 = pkin(5) + r_i_i_C(1);
t206 = sin(qJ(5));
t208 = sin(qJ(1));
t211 = cos(qJ(1));
t207 = sin(qJ(4));
t223 = qJD(1) * t207 + qJD(5);
t219 = t223 * t211;
t209 = cos(qJ(5));
t231 = qJD(5) * t207;
t224 = qJD(1) + t231;
t221 = t224 * t209;
t210 = cos(qJ(4));
t233 = qJD(4) * t210;
t201 = (t208 * t233 + t219) * t206 + t208 * t221;
t237 = t209 * pkin(5);
t203 = pkin(4) + t237;
t229 = t210 * qJD(6);
t240 = pkin(5) * t206;
t236 = r_i_i_C(3) + qJ(6) + pkin(8);
t244 = t236 * t207;
t250 = (t203 * t210 + t244) * qJD(4) - (pkin(7) - qJ(2) + t240) * qJD(1) - t231 * t240 + qJD(3) - t229;
t239 = r_i_i_C(2) * t206;
t218 = r_i_i_C(1) * t209 + t203 - t239;
t247 = t218 * t210 + t244;
t243 = r_i_i_C(2) * t209 + t241 * t206;
t234 = qJD(1) * t208;
t204 = qJD(1) * t211;
t232 = qJD(4) * t211;
t230 = qJD(5) * t210;
t225 = t236 * t210;
t220 = t224 * t206;
t217 = t243 * t207;
t215 = t223 * t208 - t210 * t232;
t213 = -qJD(5) * t237 + qJD(2) + (-t203 * t207 - pkin(1) - qJ(3) + t225) * qJD(1);
t199 = t215 * t206 - t211 * t221;
t212 = qJD(6) * t207 - t243 * t230 + (-t218 * t207 + t225) * qJD(4);
t202 = -t209 * t219 + (-t209 * t233 + t220) * t208;
t200 = t215 * t209 + t211 * t220;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t250 * t208 + t213 * t211, t204, -t234, t212 * t211 - t234 * t247, t200 * r_i_i_C(2) + t241 * t199, t207 * t232 + t210 * t234; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t213 * t208 + t250 * t211, t234, t204, t247 * t204 + t212 * t208, t202 * r_i_i_C(2) - t241 * t201, t208 * qJD(4) * t207 - t210 * t204; 0, 0, 0, -qJD(4) * t247 + qJD(5) * t217 + t229 (-t241 * t209 + t239) * t230 + qJD(4) * t217, t233;];
JaD_transl  = t1;
