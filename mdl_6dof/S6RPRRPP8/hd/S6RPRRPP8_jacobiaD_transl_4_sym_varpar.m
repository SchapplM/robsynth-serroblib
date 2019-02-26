% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:00
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPP8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPP8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP8_jacobiaD_transl_4_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:00:26
% EndTime: 2019-02-26 21:00:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (87->40), mult. (278->71), div. (0->0), fcn. (217->6), ass. (0->32)
t202 = sin(qJ(3));
t205 = cos(qJ(3));
t224 = pkin(8) + r_i_i_C(3);
t216 = t224 * t205;
t230 = -pkin(3) * t202 - qJ(2) + t216;
t201 = sin(qJ(4));
t204 = cos(qJ(4));
t210 = r_i_i_C(1) * t204 - r_i_i_C(2) * t201 + pkin(3);
t217 = t224 * t202;
t229 = t210 * t205 + t217;
t206 = cos(qJ(1));
t212 = qJD(1) * t202 + qJD(4);
t227 = t212 * t206;
t203 = sin(qJ(1));
t219 = qJD(3) * t206;
t226 = t212 * t203 - t205 * t219;
t225 = -pkin(1) - pkin(7);
t223 = qJD(1) * t203;
t222 = qJD(1) * t206;
t221 = qJD(3) * t202;
t220 = qJD(3) * t205;
t218 = qJD(4) * t205;
t213 = -qJD(4) * t202 - qJD(1);
t211 = r_i_i_C(1) * t201 + r_i_i_C(2) * t204;
t209 = t213 * t206;
t208 = t211 * qJD(4);
t207 = qJD(2) + (pkin(3) * t205 + t217) * qJD(3);
t200 = t204 * t227 + (t213 * t201 + t204 * t220) * t203;
t199 = t213 * t204 * t203 + (-t203 * t220 - t227) * t201;
t198 = t201 * t209 - t226 * t204;
t197 = t226 * t201 + t204 * t209;
t1 = [t198 * r_i_i_C(1) + t197 * r_i_i_C(2) + t207 * t206 + (t230 * t203 + t225 * t206) * qJD(1), t222, t229 * t222 + (-t211 * t218 + (-t210 * t202 + t216) * qJD(3)) * t203, r_i_i_C(1) * t199 - r_i_i_C(2) * t200, 0, 0; t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t207 * t203 + (t225 * t203 - t230 * t206) * qJD(1), t223 (t210 * t219 + t224 * t223) * t202 + (t210 * t223 + (-t224 * qJD(3) + t208) * t206) * t205, -r_i_i_C(1) * t197 + r_i_i_C(2) * t198, 0, 0; 0, 0, -t229 * qJD(3) + t202 * t208 (t201 * t218 + t204 * t221) * r_i_i_C(2) + (t201 * t221 - t204 * t218) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
