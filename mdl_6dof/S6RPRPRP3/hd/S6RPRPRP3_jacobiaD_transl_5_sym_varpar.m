% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_transl_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:42
% EndTime: 2019-02-26 20:44:42
% DurationCPUTime: 0.17s
% Computational Cost: add. (262->46), mult. (307->74), div. (0->0), fcn. (244->10), ass. (0->36)
t216 = sin(qJ(3));
t207 = cos(pkin(10)) * pkin(4) + pkin(3);
t212 = pkin(10) + qJ(5);
t208 = sin(t212);
t210 = cos(t212);
t222 = r_i_i_C(1) * t210 - r_i_i_C(2) * t208 + t207;
t217 = cos(qJ(3));
t239 = r_i_i_C(3) + pkin(8) + qJ(4);
t240 = t239 * t217;
t219 = -t222 * t216 + t240;
t244 = qJD(1) * t219;
t231 = t216 * qJD(4);
t243 = (-t207 * t216 + t240) * qJD(3) + t231;
t213 = qJ(1) + pkin(9);
t211 = cos(t213);
t237 = t210 * t211;
t236 = qJD(1) * t216;
t235 = qJD(3) * t216;
t234 = qJD(3) * t217;
t233 = qJD(5) * t216;
t232 = qJD(5) * t217;
t230 = pkin(4) * sin(pkin(10)) + pkin(7);
t229 = t239 * t216;
t226 = -qJD(1) + t232;
t225 = qJD(1) * t217 - qJD(5);
t224 = r_i_i_C(1) * t208 + r_i_i_C(2) * t210;
t223 = t226 * t208;
t221 = -t207 * t217 - pkin(2) - t229;
t209 = sin(t213);
t220 = t225 * t209 + t211 * t235;
t218 = qJD(4) * t217 + t224 * t233 + (-t222 * t217 - t229) * qJD(3);
t206 = -t225 * t237 + (t210 * t235 + t223) * t209;
t205 = t226 * t210 * t209 + (-t209 * t235 + t225 * t211) * t208;
t204 = t220 * t210 + t211 * t223;
t203 = t220 * t208 - t226 * t237;
t1 = [t206 * r_i_i_C(1) + t205 * r_i_i_C(2) - t243 * t209 + (-cos(qJ(1)) * pkin(1) - t230 * t209 + t221 * t211) * qJD(1), 0, -t209 * t244 + t218 * t211, -t209 * t236 + t211 * t234, t203 * r_i_i_C(1) + t204 * r_i_i_C(2), 0; -t204 * r_i_i_C(1) + t203 * r_i_i_C(2) + t243 * t211 + (-sin(qJ(1)) * pkin(1) + t230 * t211 + t221 * t209) * qJD(1), 0, t218 * t209 + t211 * t244, t209 * t234 + t211 * t236, -t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0; 0, 0, t219 * qJD(3) - t224 * t232 + t231, t235 (t208 * t233 - t210 * t234) * r_i_i_C(2) + (-t208 * t234 - t210 * t233) * r_i_i_C(1), 0;];
JaD_transl  = t1;
