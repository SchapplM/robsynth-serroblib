% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:58
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP5_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:58:39
% EndTime: 2019-02-26 20:58:39
% DurationCPUTime: 0.14s
% Computational Cost: add. (160->43), mult. (276->76), div. (0->0), fcn. (217->7), ass. (0->33)
t209 = pkin(9) + qJ(3);
t207 = sin(t209);
t208 = cos(t209);
t234 = pkin(8) + r_i_i_C(3);
t224 = t234 * t208;
t235 = -pkin(3) * t207 + t224;
t213 = cos(qJ(4));
t214 = cos(qJ(1));
t232 = t213 * t214;
t212 = sin(qJ(1));
t231 = qJD(1) * t212;
t230 = qJD(1) * t214;
t229 = qJD(3) * t212;
t228 = qJD(3) * t213;
t227 = qJD(3) * t214;
t226 = qJD(4) * t207;
t225 = qJD(4) * t208;
t223 = -qJD(1) + t225;
t222 = qJD(1) * t208 - qJD(4);
t211 = sin(qJ(4));
t221 = r_i_i_C(1) * t211 + r_i_i_C(2) * t213;
t220 = r_i_i_C(1) * t213 - r_i_i_C(2) * t211 + pkin(3);
t219 = t223 * t211;
t218 = -pkin(3) * t208 - t234 * t207 - cos(pkin(9)) * pkin(2) - pkin(1);
t217 = qJD(3) * t220;
t216 = t207 * t227 + t222 * t212;
t215 = -t234 * qJD(3) + t221 * qJD(4);
t210 = -pkin(7) - qJ(2);
t205 = -t222 * t232 + (t207 * t228 + t219) * t212;
t204 = t223 * t213 * t212 + (-t207 * t229 + t222 * t214) * t211;
t203 = t216 * t213 + t214 * t219;
t202 = t216 * t211 - t223 * t232;
t1 = [t205 * r_i_i_C(1) + t204 * r_i_i_C(2) + t214 * qJD(2) - t235 * t229 + (t210 * t212 + t218 * t214) * qJD(1), t230 (-t214 * t217 - t234 * t231) * t208 + (t215 * t214 + t220 * t231) * t207, t202 * r_i_i_C(1) + t203 * r_i_i_C(2), 0, 0; -t203 * r_i_i_C(1) + t202 * r_i_i_C(2) + t212 * qJD(2) + t235 * t227 + (-t210 * t214 + t218 * t212) * qJD(1), t231 (-t212 * t217 + t234 * t230) * t208 + (t215 * t212 - t220 * t230) * t207, -t204 * r_i_i_C(1) + t205 * r_i_i_C(2), 0, 0; 0, 0, -t221 * t225 + (-t220 * t207 + t224) * qJD(3) (-t208 * t228 + t211 * t226) * r_i_i_C(2) + (-qJD(3) * t208 * t211 - t213 * t226) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
