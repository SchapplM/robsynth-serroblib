% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPR8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPR8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:48
% EndTime: 2019-02-26 21:41:48
% DurationCPUTime: 0.16s
% Computational Cost: add. (170->42), mult. (303->71), div. (0->0), fcn. (242->8), ass. (0->37)
t203 = cos(pkin(10)) * pkin(3) + pkin(2);
t209 = sin(qJ(2));
t227 = t209 * qJD(3);
t211 = cos(qJ(2));
t236 = r_i_i_C(3) + pkin(8) + qJ(3);
t238 = t236 * t211;
t242 = (-t203 * t209 + t238) * qJD(2) + t227;
t206 = pkin(10) + qJ(4);
t204 = sin(t206);
t205 = cos(t206);
t217 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t203;
t214 = -t217 * t209 + t238;
t212 = cos(qJ(1));
t228 = qJD(4) * t211;
t221 = -qJD(1) + t228;
t240 = t212 * t221;
t210 = sin(qJ(1));
t220 = qJD(1) * t211 - qJD(4);
t232 = qJD(2) * t209;
t237 = -t210 * t232 + t220 * t212;
t234 = qJD(1) * t210;
t233 = qJD(1) * t212;
t231 = qJD(2) * t211;
t230 = qJD(2) * t212;
t229 = qJD(4) * t209;
t226 = pkin(3) * sin(pkin(10)) + pkin(7);
t224 = t236 * t209;
t219 = r_i_i_C(1) * t204 + r_i_i_C(2) * t205;
t218 = t221 * t210;
t216 = -t203 * t211 - pkin(1) - t224;
t215 = t209 * t230 + t220 * t210;
t213 = qJD(3) * t211 + t219 * t229 + (-t217 * t211 - t224) * qJD(2);
t202 = t204 * t218 - t237 * t205;
t201 = t237 * t204 + t205 * t218;
t200 = t204 * t240 + t215 * t205;
t199 = t215 * t204 - t205 * t240;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) - t242 * t210 + (-t226 * t210 + t216 * t212) * qJD(1), t213 * t212 - t214 * t234, -t209 * t234 + t211 * t230, t199 * r_i_i_C(1) + t200 * r_i_i_C(2), 0, 0; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t242 * t212 + (t216 * t210 + t226 * t212) * qJD(1), t213 * t210 + t214 * t233, t209 * t233 + t210 * t231, -t201 * r_i_i_C(1) + t202 * r_i_i_C(2), 0, 0; 0, t214 * qJD(2) - t219 * t228 + t227, t232 (t204 * t229 - t205 * t231) * r_i_i_C(2) + (-t204 * t231 - t205 * t229) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
