% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP6
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
% Datum: 2019-02-26 20:46
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RPRPRP6_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP6_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RPRPRP6_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP6_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:46:28
% EndTime: 2019-02-26 20:46:28
% DurationCPUTime: 0.18s
% Computational Cost: add. (197->47), mult. (324->78), div. (0->0), fcn. (254->7), ass. (0->35)
t205 = pkin(9) + qJ(3);
t203 = sin(t205);
t204 = cos(t205);
t223 = pkin(3) + pkin(8) + r_i_i_C(3);
t220 = t223 * t203;
t236 = (-qJ(4) * t204 + t220) * qJD(3) - t203 * qJD(4);
t209 = cos(qJ(5));
t219 = qJD(5) * t203 + qJD(1);
t234 = t209 * t219;
t207 = sin(qJ(5));
t233 = t219 * t207;
t232 = pkin(4) + pkin(7) + qJ(2);
t208 = sin(qJ(1));
t230 = qJD(1) * t208;
t210 = cos(qJ(1));
t229 = qJD(1) * t210;
t228 = qJD(3) * t203;
t227 = qJD(3) * t204;
t226 = qJD(3) * t209;
t225 = qJD(5) * t204;
t222 = t208 * t227;
t221 = t210 * t227;
t218 = -qJD(1) * t203 - qJD(5);
t217 = t218 * t210;
t216 = r_i_i_C(1) * t207 + r_i_i_C(2) * t209 + qJ(4);
t215 = qJD(3) * t216;
t214 = -qJ(4) * t203 - t223 * t204 - cos(pkin(9)) * pkin(2) - pkin(1);
t213 = qJD(4) + (r_i_i_C(1) * t209 - r_i_i_C(2) * t207) * qJD(5);
t212 = t218 * t208 + t221;
t211 = -t223 * qJD(3) + t213;
t201 = t212 * t207 + t210 * t234;
t200 = t212 * t209 - t210 * t233;
t199 = -t208 * t234 + (t217 - t222) * t207;
t198 = t209 * t217 + (-t204 * t226 + t233) * t208;
t1 = [t199 * r_i_i_C(1) + t198 * r_i_i_C(2) + t210 * qJD(2) + t236 * t208 + (-t232 * t208 + t214 * t210) * qJD(1), t229 (-t210 * t215 + t223 * t230) * t203 + (t211 * t210 - t216 * t230) * t204, -t203 * t230 + t221, t200 * r_i_i_C(1) - t201 * r_i_i_C(2), 0; t201 * r_i_i_C(1) + t200 * r_i_i_C(2) + t208 * qJD(2) - t236 * t210 + (t214 * t208 + t232 * t210) * qJD(1), t230 (-t208 * t215 - t223 * t229) * t203 + (t211 * t208 + t216 * t229) * t204, t203 * t229 + t222, -t198 * r_i_i_C(1) + t199 * r_i_i_C(2), 0; 0, 0, t213 * t203 + (t216 * t204 - t220) * qJD(3), t228 (-t207 * t228 + t209 * t225) * r_i_i_C(2) + (t203 * t226 + t207 * t225) * r_i_i_C(1), 0;];
JaD_transl  = t1;
