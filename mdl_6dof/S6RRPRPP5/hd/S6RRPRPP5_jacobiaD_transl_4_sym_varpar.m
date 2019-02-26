% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_transl_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:08
% EndTime: 2019-02-26 21:37:08
% DurationCPUTime: 0.17s
% Computational Cost: add. (101->43), mult. (318->73), div. (0->0), fcn. (248->6), ass. (0->34)
t199 = sin(qJ(2));
t202 = cos(qJ(2));
t216 = pkin(2) + pkin(8) + r_i_i_C(3);
t213 = t216 * t199;
t229 = (-qJ(3) * t202 + t213) * qJD(2) - t199 * qJD(3);
t201 = cos(qJ(4));
t211 = qJD(4) * t199 + qJD(1);
t227 = t201 * t211;
t198 = sin(qJ(4));
t226 = t211 * t198;
t225 = pkin(3) + pkin(7);
t200 = sin(qJ(1));
t223 = qJD(1) * t200;
t203 = cos(qJ(1));
t222 = qJD(1) * t203;
t221 = qJD(2) * t199;
t220 = qJD(2) * t202;
t219 = qJD(2) * t203;
t218 = qJD(4) * t202;
t215 = t200 * t220;
t214 = t202 * t219;
t212 = t216 * t202;
t210 = -qJD(1) * t199 - qJD(4);
t209 = t210 * t203;
t208 = r_i_i_C(1) * t198 + r_i_i_C(2) * t201 + qJ(3);
t207 = -qJ(3) * t199 - pkin(1) - t212;
t206 = qJD(3) + (r_i_i_C(1) * t201 - r_i_i_C(2) * t198) * qJD(4);
t205 = t210 * t200 + t214;
t204 = t208 * t202 - t213;
t197 = t205 * t198 + t203 * t227;
t196 = t205 * t201 - t203 * t226;
t195 = -t200 * t227 + (t209 - t215) * t198;
t194 = t201 * t209 + (-t201 * t220 + t226) * t200;
t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t229 * t200 + (-t225 * t200 + t207 * t203) * qJD(1) (-t208 * t219 + t216 * t223) * t199 + (-t208 * t223 + (-t216 * qJD(2) + t206) * t203) * t202, -t199 * t223 + t214, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0, 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t229 * t203 + (t207 * t200 + t225 * t203) * qJD(1), t204 * t222 + (t206 * t202 + (-t208 * t199 - t212) * qJD(2)) * t200, t199 * t222 + t215, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; 0, t204 * qJD(2) + t206 * t199, t221 (-t198 * t221 + t201 * t218) * r_i_i_C(2) + (t198 * t218 + t201 * t221) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
