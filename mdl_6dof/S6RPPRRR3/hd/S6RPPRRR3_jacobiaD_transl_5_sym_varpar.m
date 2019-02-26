% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:35:57
% EndTime: 2019-02-26 20:35:58
% DurationCPUTime: 0.15s
% Computational Cost: add. (175->40), mult. (282->67), div. (0->0), fcn. (219->8), ass. (0->33)
t202 = sin(qJ(4));
t201 = sin(qJ(5));
t203 = cos(qJ(5));
t210 = r_i_i_C(1) * t203 - r_i_i_C(2) * t201 + pkin(4);
t204 = cos(qJ(4));
t221 = pkin(8) + r_i_i_C(3);
t223 = t221 * t204;
t231 = (-t210 * t202 + t223) * qJD(4);
t216 = t221 * t202;
t230 = t210 * t204 + t216;
t229 = -pkin(4) * t202 - qJ(3) + t223;
t212 = qJD(1) * t202 + qJD(5);
t227 = t201 * t212;
t226 = t203 * t212;
t222 = -pkin(2) - pkin(7);
t220 = qJD(4) * t202;
t219 = qJD(4) * t204;
t218 = qJD(5) * t204;
t213 = -qJD(5) * t202 - qJD(1);
t211 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203;
t209 = qJD(5) * t211;
t208 = qJD(3) + (pkin(4) * t204 + t216) * qJD(4);
t207 = -t201 * t219 + t213 * t203;
t206 = t213 * t201 + t203 * t219;
t205 = qJD(1) * t230;
t200 = qJ(1) + pkin(10);
t199 = cos(t200);
t198 = sin(t200);
t197 = t206 * t198 + t199 * t226;
t196 = t207 * t198 - t199 * t227;
t195 = -t198 * t226 + t206 * t199;
t194 = t198 * t227 + t207 * t199;
t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t208 * t199 + (-cos(qJ(1)) * pkin(1) + t222 * t199 + t229 * t198) * qJD(1), 0, qJD(1) * t199, t199 * t205 + (-t211 * t218 + t231) * t198, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t208 * t198 + (-sin(qJ(1)) * pkin(1) + t222 * t198 - t229 * t199) * qJD(1), 0, qJD(1) * t198, t198 * t205 + (t204 * t209 - t231) * t199, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, 0, -t230 * qJD(4) + t202 * t209 (t201 * t218 + t203 * t220) * r_i_i_C(2) + (t201 * t220 - t203 * t218) * r_i_i_C(1), 0;];
JaD_transl  = t1;
