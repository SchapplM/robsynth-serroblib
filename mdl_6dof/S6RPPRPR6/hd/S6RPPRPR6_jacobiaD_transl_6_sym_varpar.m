% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:28
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPPRPR6_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR6_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPPRPR6_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRPR6_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:28:36
% EndTime: 2019-02-26 20:28:36
% DurationCPUTime: 0.21s
% Computational Cost: add. (110->48), mult. (334->77), div. (0->0), fcn. (260->6), ass. (0->33)
t198 = sin(qJ(4));
t201 = cos(qJ(4));
t216 = pkin(4) + pkin(8) + r_i_i_C(3);
t212 = t216 * t201;
t226 = (qJ(5) * t198 + t212) * qJD(4) - t201 * qJD(5) + qJD(3);
t202 = cos(qJ(1));
t220 = qJD(1) * t201;
t210 = qJD(6) + t220;
t224 = t210 * t202;
t199 = sin(qJ(1));
t219 = qJD(4) * t198;
t213 = t202 * t219;
t223 = t210 * t199 + t213;
t221 = qJD(1) * t199;
t196 = qJD(1) * t202;
t218 = qJD(4) * t201;
t217 = qJD(6) * t198;
t215 = pkin(5) + pkin(7) - qJ(2);
t214 = t199 * t219;
t211 = qJD(6) * t201 + qJD(1);
t208 = t211 * t202;
t197 = sin(qJ(6));
t200 = cos(qJ(6));
t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200 + qJ(5);
t206 = qJD(4) * t207;
t205 = qJD(5) + (r_i_i_C(1) * t200 - r_i_i_C(2) * t197) * qJD(6);
t204 = qJ(5) * t201 - t216 * t198 - pkin(1) - qJ(3);
t203 = -t216 * qJD(4) + t205;
t195 = -t223 * t197 + t200 * t208;
t194 = t197 * t208 + t223 * t200;
t193 = t211 * t200 * t199 + (-t214 + t224) * t197;
t192 = -t200 * t224 + (t211 * t197 + t200 * t219) * t199;
t1 = [t193 * r_i_i_C(1) - t192 * r_i_i_C(2) + t202 * qJD(2) - t226 * t199 + (t215 * t199 + t204 * t202) * qJD(1), t196, -t221 (t202 * t206 - t216 * t221) * t201 + (t203 * t202 - t207 * t221) * t198, t199 * t220 + t213, t194 * r_i_i_C(1) + t195 * r_i_i_C(2); -t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t199 * qJD(2) + t226 * t202 + (t204 * t199 - t215 * t202) * qJD(1), t221, t196 (t216 * t196 + t199 * t206) * t201 + (t207 * t196 + t203 * t199) * t198, -t201 * t196 + t214, t192 * r_i_i_C(1) + t193 * r_i_i_C(2); 0, 0, 0, t205 * t201 + (-t207 * t198 - t212) * qJD(4), t218 (-t197 * t218 - t200 * t217) * r_i_i_C(2) + (-t197 * t217 + t200 * t218) * r_i_i_C(1);];
JaD_transl  = t1;
