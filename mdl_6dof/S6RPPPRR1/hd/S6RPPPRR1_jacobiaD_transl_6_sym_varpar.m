% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:22
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPPRR1_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPPRR1_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:22:35
% EndTime: 2019-02-26 20:22:36
% DurationCPUTime: 0.15s
% Computational Cost: add. (186->43), mult. (290->66), div. (0->0), fcn. (225->8), ass. (0->34)
t206 = cos(qJ(5));
t204 = sin(qJ(5));
t227 = pkin(8) + r_i_i_C(3);
t229 = t227 * t204;
t234 = (pkin(5) * t206 + t229) * qJD(5) + qJD(4);
t203 = sin(qJ(6));
t205 = cos(qJ(6));
t212 = r_i_i_C(1) * t205 - r_i_i_C(2) * t203 + pkin(5);
t232 = t212 * t206 + t229;
t214 = qJD(1) * t204 + qJD(6);
t231 = t205 * t214;
t221 = qJD(6) * t204;
t215 = qJD(1) + t221;
t222 = qJD(5) * t206;
t228 = t203 * t222 + t215 * t205;
t225 = pkin(7) - qJ(3);
t202 = qJ(1) + pkin(9);
t200 = sin(t202);
t224 = qJD(1) * t200;
t201 = cos(t202);
t199 = qJD(1) * t201;
t223 = qJD(5) * t204;
t220 = qJD(6) * t206;
t217 = t227 * t206;
t213 = r_i_i_C(1) * t203 + r_i_i_C(2) * t205;
t211 = t214 * t203;
t210 = -pkin(5) * t204 - pkin(2) - qJ(4) + t217;
t209 = t215 * t203 - t205 * t222;
t207 = -t213 * t220 + (-t212 * t204 + t217) * qJD(5);
t198 = t209 * t200 - t201 * t231;
t197 = t228 * t200 + t201 * t211;
t196 = t200 * t231 + t209 * t201;
t195 = t200 * t211 - t228 * t201;
t1 = [t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t201 * qJD(3) - t234 * t200 + (-cos(qJ(1)) * pkin(1) + t225 * t200 + t210 * t201) * qJD(1), 0, t199, -t224, t207 * t201 - t224 * t232, t195 * r_i_i_C(1) + t196 * r_i_i_C(2); -t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t200 * qJD(3) + t234 * t201 + (-sin(qJ(1)) * pkin(1) - t225 * t201 + t210 * t200) * qJD(1), 0, t224, t199, t232 * t199 + t207 * t200, -t197 * r_i_i_C(1) + t198 * r_i_i_C(2); 0, 0, 0, 0, -qJD(5) * t232 + t213 * t221 (t203 * t220 + t205 * t223) * r_i_i_C(2) + (t203 * t223 - t205 * t220) * r_i_i_C(1);];
JaD_transl  = t1;
