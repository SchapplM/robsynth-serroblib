% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 3 (0=Basis) von
% S6RRRRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:48
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRRRRR4_jacobiaD_transl_3_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_transl_3_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRR4_jacobiaD_transl_3_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRRRRR4_jacobiaD_transl_3_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRR4_jacobiaD_transl_3_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_3_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:48:48
% EndTime: 2019-02-26 22:48:48
% DurationCPUTime: 0.14s
% Computational Cost: add. (83->38), mult. (270->73), div. (0->0), fcn. (211->6), ass. (0->31)
t199 = sin(qJ(2));
t202 = cos(qJ(2));
t223 = pkin(8) + r_i_i_C(3);
t213 = t223 * t202;
t224 = -pkin(2) * t199 + t213;
t201 = cos(qJ(3));
t203 = cos(qJ(1));
t221 = t201 * t203;
t200 = sin(qJ(1));
t220 = qJD(1) * t200;
t219 = qJD(1) * t203;
t218 = qJD(2) * t200;
t217 = qJD(2) * t202;
t216 = qJD(2) * t203;
t215 = qJD(3) * t199;
t214 = qJD(3) * t202;
t212 = -qJD(1) + t214;
t211 = qJD(1) * t202 - qJD(3);
t198 = sin(qJ(3));
t210 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201;
t209 = r_i_i_C(1) * t201 - r_i_i_C(2) * t198 + pkin(2);
t208 = t212 * t198;
t207 = -pkin(2) * t202 - t223 * t199 - pkin(1);
t206 = qJD(2) * t209;
t205 = t199 * t216 + t211 * t200;
t204 = -t223 * qJD(2) + t210 * qJD(3);
t197 = -t211 * t221 + (qJD(2) * t199 * t201 + t208) * t200;
t196 = t212 * t201 * t200 + (-t199 * t218 + t211 * t203) * t198;
t195 = t205 * t201 + t203 * t208;
t194 = t205 * t198 - t212 * t221;
t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t224 * t218 + (-pkin(7) * t200 + t207 * t203) * qJD(1) (-t203 * t206 - t223 * t220) * t202 + (t204 * t203 + t209 * t220) * t199, t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0, 0; -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t224 * t216 + (pkin(7) * t203 + t207 * t200) * qJD(1) (-t200 * t206 + t223 * t219) * t202 + (t204 * t200 - t209 * t219) * t199, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0, 0; 0, -t210 * t214 + (-t209 * t199 + t213) * qJD(2) (t198 * t215 - t201 * t217) * r_i_i_C(2) + (-t198 * t217 - t201 * t215) * r_i_i_C(1), 0, 0, 0;];
JaD_transl  = t1;
