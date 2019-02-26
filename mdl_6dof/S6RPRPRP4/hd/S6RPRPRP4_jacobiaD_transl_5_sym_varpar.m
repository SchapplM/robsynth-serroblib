% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:45
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:45:18
% EndTime: 2019-02-26 20:45:19
% DurationCPUTime: 0.23s
% Computational Cost: add. (201->42), mult. (322->69), div. (0->0), fcn. (250->8), ass. (0->33)
t204 = cos(qJ(3));
t201 = sin(qJ(5));
t203 = cos(qJ(5));
t209 = r_i_i_C(1) * t201 + r_i_i_C(2) * t203 + qJ(4);
t202 = sin(qJ(3));
t219 = pkin(3) + pkin(8) + r_i_i_C(3);
t229 = t219 * t202;
t231 = -t209 * t204 + t229;
t234 = qJD(1) * t231;
t233 = (-qJ(4) * t204 + t229) * qJD(3) - t202 * qJD(4);
t213 = qJD(5) * t202 + qJD(1);
t222 = qJD(3) * t204;
t228 = t213 * t201 - t203 * t222;
t227 = t201 * t222 + t213 * t203;
t226 = pkin(4) + pkin(7);
t224 = qJD(1) * t202;
t223 = qJD(3) * t202;
t221 = qJD(5) * t204;
t215 = t219 * t204;
t212 = -qJD(5) - t224;
t211 = t212 * t201;
t210 = t212 * t203;
t208 = -qJ(4) * t202 - pkin(2) - t215;
t207 = qJD(4) + (r_i_i_C(1) * t203 - r_i_i_C(2) * t201) * qJD(5);
t205 = t207 * t204 + (-t209 * t202 - t215) * qJD(3);
t200 = qJ(1) + pkin(9);
t199 = cos(t200);
t198 = sin(t200);
t197 = t198 * t211 + t227 * t199;
t196 = t198 * t210 - t228 * t199;
t195 = -t227 * t198 + t199 * t211;
t194 = t228 * t198 + t199 * t210;
t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) + t233 * t198 + (-cos(qJ(1)) * pkin(1) - t226 * t198 + t208 * t199) * qJD(1), 0, t198 * t234 + t205 * t199, -t198 * t224 + t199 * t222, t196 * r_i_i_C(1) - t197 * r_i_i_C(2), 0; t197 * r_i_i_C(1) + t196 * r_i_i_C(2) - t233 * t199 + (-sin(qJ(1)) * pkin(1) + t226 * t199 + t208 * t198) * qJD(1), 0, t205 * t198 - t199 * t234, t198 * t222 + t199 * t224, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0; 0, 0, -qJD(3) * t231 + t207 * t202, t223 (-t201 * t223 + t203 * t221) * r_i_i_C(2) + (t201 * t221 + t203 * t223) * r_i_i_C(1), 0;];
JaD_transl  = t1;
