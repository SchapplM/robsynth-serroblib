% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR8_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR8_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:45
% EndTime: 2019-02-26 20:29:45
% DurationCPUTime: 0.19s
% Computational Cost: add. (202->49), mult. (336->81), div. (0->0), fcn. (262->7), ass. (0->36)
t208 = pkin(9) + qJ(4);
t206 = cos(t208);
t229 = pkin(4) + pkin(8) + r_i_i_C(3);
t239 = t229 * t206;
t205 = sin(t208);
t238 = t229 * t205 + pkin(3) * sin(pkin(9)) - qJ(5) * t206 + qJ(2);
t211 = sin(qJ(6));
t213 = cos(qJ(6));
t217 = qJD(5) + (r_i_i_C(1) * t213 - r_i_i_C(2) * t211) * qJD(6);
t237 = -t229 * qJD(4) + t217;
t214 = cos(qJ(1));
t236 = t213 * t214;
t212 = sin(qJ(1));
t235 = qJD(1) * t212;
t207 = qJD(1) * t214;
t234 = qJD(4) * t206;
t233 = qJD(4) * t212;
t232 = qJD(4) * t213;
t231 = qJD(4) * t214;
t230 = qJD(6) * t205;
t228 = -pkin(1) - pkin(5) - pkin(7) - qJ(3);
t227 = t205 * t233;
t226 = t205 * t231;
t225 = qJD(1) * t229;
t224 = qJD(6) * t206 + qJD(1);
t223 = qJD(1) * t206 + qJD(6);
t221 = t224 * t211;
t220 = r_i_i_C(1) * t211 + r_i_i_C(2) * t213 + qJ(5);
t218 = qJD(1) * t220;
t216 = t223 * t212 + t226;
t215 = -t206 * qJD(5) + qJD(2) + (qJ(5) * t205 + t239) * qJD(4);
t204 = t216 * t211 - t224 * t236;
t203 = t216 * t213 + t214 * t221;
t202 = t224 * t213 * t212 + (t223 * t214 - t227) * t211;
t201 = -t223 * t236 + (t205 * t232 + t221) * t212;
t1 = [t204 * r_i_i_C(1) + t203 * r_i_i_C(2) - t212 * qJD(3) + t215 * t214 + (-t238 * t212 + t228 * t214) * qJD(1), t207, -t235 (t214 * t225 + t220 * t233) * t206 + (t237 * t212 + t214 * t218) * t205, -t206 * t207 + t227, t201 * r_i_i_C(1) + t202 * r_i_i_C(2); -t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t214 * qJD(3) + t215 * t212 + (t228 * t212 + t238 * t214) * qJD(1), t235, t207 (t212 * t225 - t220 * t231) * t206 + (t212 * t218 - t237 * t214) * t205, -t206 * t235 - t226, -t203 * r_i_i_C(1) + t204 * r_i_i_C(2); 0, 0, 0, t217 * t206 + (-t220 * t205 - t239) * qJD(4), t234 (-t211 * t234 - t213 * t230) * r_i_i_C(2) + (t206 * t232 - t211 * t230) * r_i_i_C(1);];
JaD_transl  = t1;
