% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:26:48
% EndTime: 2019-02-26 20:26:48
% DurationCPUTime: 0.19s
% Computational Cost: add. (272->46), mult. (312->71), div. (0->0), fcn. (240->10), ass. (0->36)
t220 = qJ(4) + pkin(10);
t216 = sin(t220);
t218 = cos(t220);
t247 = pkin(8) + r_i_i_C(3);
t234 = t247 * t218 - sin(qJ(4)) * pkin(4);
t223 = sin(qJ(6));
t225 = cos(qJ(6));
t235 = r_i_i_C(1) * t225 - r_i_i_C(2) * t223 + pkin(5);
t256 = (-t235 * t216 + t234) * qJD(4);
t255 = -pkin(5) * t216 - qJ(3) + t234;
t232 = t247 * t216 + cos(qJ(4)) * pkin(4);
t253 = t235 * t218 + t232;
t237 = qJD(1) * t216 + qJD(6);
t252 = t223 * t237;
t251 = t225 * t237;
t243 = -pkin(2) - qJ(5) - pkin(7);
t221 = qJ(1) + pkin(9);
t217 = sin(t221);
t242 = qJD(1) * t217;
t241 = qJD(4) * t223;
t240 = qJD(4) * t225;
t239 = qJD(6) * t218;
t238 = -qJD(6) * t216 - qJD(1);
t236 = r_i_i_C(1) * t223 + r_i_i_C(2) * t225;
t231 = qJD(6) * t236;
t230 = -t218 * t241 + t238 * t225;
t229 = t218 * t240 + t238 * t223;
t228 = qJD(3) + (pkin(5) * t218 + t232) * qJD(4);
t227 = qJD(1) * t253;
t219 = cos(t221);
t215 = qJD(1) * t219;
t214 = t229 * t217 + t219 * t251;
t213 = t230 * t217 - t219 * t252;
t212 = -t217 * t251 + t229 * t219;
t211 = t217 * t252 + t230 * t219;
t1 = [t212 * r_i_i_C(1) + t211 * r_i_i_C(2) - t217 * qJD(5) + t228 * t219 + (-cos(qJ(1)) * pkin(1) + t243 * t219 + t255 * t217) * qJD(1), 0, t215, t219 * t227 + (-t236 * t239 + t256) * t217, -t242, t213 * r_i_i_C(1) - t214 * r_i_i_C(2); t214 * r_i_i_C(1) + t213 * r_i_i_C(2) + t219 * qJD(5) + t228 * t217 + (-sin(qJ(1)) * pkin(1) + t243 * t217 - t255 * t219) * qJD(1), 0, t242, t217 * t227 + (t218 * t231 - t256) * t219, t215, -t211 * r_i_i_C(1) + t212 * r_i_i_C(2); 0, 0, 0, -t253 * qJD(4) + t216 * t231, 0 (t216 * t240 + t223 * t239) * r_i_i_C(2) + (t216 * t241 - t225 * t239) * r_i_i_C(1);];
JaD_transl  = t1;
