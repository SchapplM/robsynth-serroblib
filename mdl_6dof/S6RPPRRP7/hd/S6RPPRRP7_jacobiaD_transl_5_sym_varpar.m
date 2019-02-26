% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRP7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:33
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRP7_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRP7_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:33:37
% EndTime: 2019-02-26 20:33:37
% DurationCPUTime: 0.16s
% Computational Cost: add. (165->45), mult. (288->78), div. (0->0), fcn. (225->7), ass. (0->33)
t211 = sin(qJ(5));
t213 = cos(qJ(5));
t217 = (r_i_i_C(1) * t211 + r_i_i_C(2) * t213) * qJD(5);
t232 = pkin(8) + r_i_i_C(3);
t238 = t232 * qJD(4) - t217;
t214 = cos(qJ(1));
t208 = pkin(9) + qJ(4);
t205 = sin(t208);
t222 = qJD(1) * t205 + qJD(5);
t237 = t222 * t214;
t236 = t232 * t205;
t206 = cos(t208);
t234 = t232 * t206 - pkin(3) * sin(pkin(9)) - pkin(4) * t205 - qJ(2);
t212 = sin(qJ(1));
t227 = qJD(4) * t214;
t233 = -t206 * t227 + t222 * t212;
t231 = -pkin(1) - pkin(7) - qJ(3);
t230 = qJD(1) * t212;
t229 = qJD(4) * t212;
t228 = qJD(4) * t213;
t226 = qJD(5) * t206;
t225 = qJD(1) * t232;
t223 = -qJD(5) * t205 - qJD(1);
t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(4);
t219 = t223 * t214;
t216 = qJD(1) * t220;
t215 = qJD(2) + (pkin(4) * t206 + t236) * qJD(4);
t207 = qJD(1) * t214;
t204 = t213 * t237 + (t206 * t228 + t223 * t211) * t212;
t203 = t223 * t213 * t212 + (-t206 * t229 - t237) * t211;
t202 = t211 * t219 - t233 * t213;
t201 = t233 * t211 + t213 * t219;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t212 * qJD(3) + t215 * t214 + (t234 * t212 + t231 * t214) * qJD(1), t207, -t230 (t214 * t225 - t220 * t229) * t205 + (t238 * t212 + t214 * t216) * t206, t203 * r_i_i_C(1) - t204 * r_i_i_C(2), 0; t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t214 * qJD(3) + t215 * t212 + (t231 * t212 - t234 * t214) * qJD(1), t230, t207 (t212 * t225 + t220 * t227) * t205 + (t212 * t216 - t238 * t214) * t206, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0; 0, 0, 0, t205 * t217 + (-t220 * t206 - t236) * qJD(4) (t205 * t228 + t211 * t226) * r_i_i_C(2) + (qJD(4) * t205 * t211 - t213 * t226) * r_i_i_C(1), 0;];
JaD_transl  = t1;
