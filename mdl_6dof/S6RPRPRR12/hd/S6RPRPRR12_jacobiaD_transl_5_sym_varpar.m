% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:55
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRR12_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRR12_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:55:13
% EndTime: 2019-02-26 20:55:13
% DurationCPUTime: 0.19s
% Computational Cost: add. (105->44), mult. (326->75), div. (0->0), fcn. (254->6), ass. (0->34)
t201 = cos(qJ(3));
t216 = pkin(3) + pkin(8) + r_i_i_C(3);
t227 = t216 * t201;
t198 = sin(qJ(3));
t226 = -qJ(4) * t201 + t216 * t198 + qJ(2);
t197 = sin(qJ(5));
t200 = cos(qJ(5));
t205 = qJD(4) + (r_i_i_C(1) * t200 - r_i_i_C(2) * t197) * qJD(5);
t225 = -t216 * qJD(3) + t205;
t202 = cos(qJ(1));
t224 = t200 * t202;
t199 = sin(qJ(1));
t223 = qJD(1) * t199;
t222 = qJD(1) * t201;
t221 = qJD(1) * t202;
t220 = qJD(3) * t199;
t219 = qJD(3) * t201;
t218 = qJD(3) * t202;
t217 = qJD(5) * t198;
t215 = t198 * t220;
t214 = t198 * t218;
t212 = qJD(1) * t216;
t211 = qJD(5) * t201 + qJD(1);
t210 = qJD(5) + t222;
t208 = t211 * t197;
t207 = r_i_i_C(1) * t197 + r_i_i_C(2) * t200 + qJ(4);
t206 = qJD(1) * t207;
t204 = t210 * t199 + t214;
t203 = -t201 * qJD(4) + qJD(2) + (-pkin(1) - pkin(4) - pkin(7)) * qJD(1) + (qJ(4) * t198 + t227) * qJD(3);
t196 = t204 * t197 - t211 * t224;
t195 = t204 * t200 + t202 * t208;
t194 = t211 * t200 * t199 + (t210 * t202 - t215) * t197;
t193 = -t210 * t224 + (qJD(3) * t198 * t200 + t208) * t199;
t1 = [t196 * r_i_i_C(1) + t195 * r_i_i_C(2) + t203 * t202 - t226 * t223, t221 (t202 * t212 + t207 * t220) * t201 + (t225 * t199 + t202 * t206) * t198, -t201 * t221 + t215, r_i_i_C(1) * t193 + r_i_i_C(2) * t194, 0; -t194 * r_i_i_C(1) + t193 * r_i_i_C(2) + t203 * t199 + t226 * t221, t223 (t199 * t212 - t207 * t218) * t201 + (t199 * t206 - t225 * t202) * t198, -t199 * t222 - t214, -r_i_i_C(1) * t195 + r_i_i_C(2) * t196, 0; 0, 0, t205 * t201 + (-t207 * t198 - t227) * qJD(3), t219 (-t197 * t219 - t200 * t217) * r_i_i_C(2) + (-t197 * t217 + t200 * t219) * r_i_i_C(1), 0;];
JaD_transl  = t1;
