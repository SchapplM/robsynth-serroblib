% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR7_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR7_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_jacobiaD_transl_6_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:10
% EndTime: 2019-02-26 20:29:10
% DurationCPUTime: 0.18s
% Computational Cost: add. (262->51), mult. (321->78), div. (0->0), fcn. (256->9), ass. (0->39)
t218 = pkin(9) + qJ(4);
t213 = sin(t218);
t211 = cos(pkin(10)) * pkin(5) + pkin(4);
t217 = pkin(10) + qJ(6);
t212 = sin(t217);
t214 = cos(t217);
t229 = r_i_i_C(1) * t214 - r_i_i_C(2) * t212 + t211;
t215 = cos(t218);
t248 = r_i_i_C(3) + pkin(8) + qJ(5);
t251 = t248 * t215;
t257 = -(-t229 * t213 + t251) * qJD(4) - qJD(5) * t213;
t237 = t248 * t213;
t256 = t229 * t215 + t237;
t255 = t251 - pkin(3) * sin(pkin(9)) - t211 * t213 - qJ(2);
t223 = sin(qJ(1));
t234 = qJD(1) * t213 + qJD(6);
t224 = cos(qJ(1));
t244 = qJD(4) * t224;
t250 = -t215 * t244 + t234 * t223;
t245 = qJD(4) * t223;
t249 = t215 * t245 + t234 * t224;
t247 = qJD(1) * t223;
t216 = qJD(1) * t224;
t246 = qJD(4) * t213;
t242 = qJD(6) * t215;
t241 = t215 * qJD(5);
t235 = -qJD(6) * t213 - qJD(1);
t233 = -pkin(5) * sin(pkin(10)) - pkin(1) - pkin(7) - qJ(3);
t232 = r_i_i_C(1) * t212 + r_i_i_C(2) * t214;
t231 = t235 * t223;
t230 = t235 * t224;
t227 = qJD(6) * t232;
t226 = -t241 + qJD(2) + (t211 * t215 + t237) * qJD(4);
t225 = qJD(1) * t256;
t210 = t212 * t231 + t249 * t214;
t209 = -t249 * t212 + t214 * t231;
t208 = t212 * t230 - t250 * t214;
t207 = t250 * t212 + t214 * t230;
t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) - t223 * qJD(3) + t226 * t224 + (t255 * t223 + t233 * t224) * qJD(1), t216, -t247, t224 * t225 + (-t232 * t242 - t257) * t223, t213 * t245 - t215 * t216, t209 * r_i_i_C(1) - t210 * r_i_i_C(2); t210 * r_i_i_C(1) + t209 * r_i_i_C(2) + t224 * qJD(3) + t226 * t223 + (t233 * t223 - t255 * t224) * qJD(1), t247, t216, t223 * t225 + (t215 * t227 + t257) * t224, -t213 * t244 - t215 * t247, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2); 0, 0, 0, -t256 * qJD(4) + t213 * t227 + t241, qJD(4) * t215 (t212 * t242 + t214 * t246) * r_i_i_C(2) + (t212 * t246 - t214 * t242) * r_i_i_C(1);];
JaD_transl  = t1;
