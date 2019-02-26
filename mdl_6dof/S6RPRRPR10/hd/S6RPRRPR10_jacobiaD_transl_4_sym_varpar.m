% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRPR10
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:06
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRRPR10_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRRPR10_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:06:04
% EndTime: 2019-02-26 21:06:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (87->40), mult. (278->71), div. (0->0), fcn. (217->6), ass. (0->32)
t201 = sin(qJ(3));
t204 = cos(qJ(3));
t223 = pkin(8) + r_i_i_C(3);
t215 = t223 * t204;
t229 = -pkin(3) * t201 - qJ(2) + t215;
t200 = sin(qJ(4));
t203 = cos(qJ(4));
t209 = r_i_i_C(1) * t203 - r_i_i_C(2) * t200 + pkin(3);
t216 = t223 * t201;
t228 = t209 * t204 + t216;
t205 = cos(qJ(1));
t211 = qJD(1) * t201 + qJD(4);
t226 = t211 * t205;
t202 = sin(qJ(1));
t218 = qJD(3) * t205;
t225 = t211 * t202 - t204 * t218;
t224 = -pkin(1) - pkin(7);
t222 = qJD(1) * t202;
t221 = qJD(1) * t205;
t220 = qJD(3) * t201;
t219 = qJD(3) * t204;
t217 = qJD(4) * t204;
t212 = -qJD(4) * t201 - qJD(1);
t210 = r_i_i_C(1) * t200 + r_i_i_C(2) * t203;
t208 = t212 * t205;
t207 = t210 * qJD(4);
t206 = qJD(2) + (pkin(3) * t204 + t216) * qJD(3);
t199 = t203 * t226 + (t212 * t200 + t203 * t219) * t202;
t198 = t212 * t203 * t202 + (-t202 * t219 - t226) * t200;
t197 = t200 * t208 - t225 * t203;
t196 = t225 * t200 + t203 * t208;
t1 = [t197 * r_i_i_C(1) + t196 * r_i_i_C(2) + t206 * t205 + (t229 * t202 + t224 * t205) * qJD(1), t221, t228 * t221 + (-t210 * t217 + (-t209 * t201 + t215) * qJD(3)) * t202, t198 * r_i_i_C(1) - t199 * r_i_i_C(2), 0, 0; t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t206 * t202 + (t224 * t202 - t229 * t205) * qJD(1), t222 (t209 * t218 + t223 * t222) * t201 + (t209 * t222 + (-t223 * qJD(3) + t207) * t205) * t204, -t196 * r_i_i_C(1) + t197 * r_i_i_C(2), 0, 0; 0, 0, -t228 * qJD(3) + t201 * t207 (t200 * t217 + t203 * t220) * r_i_i_C(2) + (t200 * t220 - t203 * t217) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
