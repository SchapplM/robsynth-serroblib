% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_transl_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:09
% EndTime: 2019-02-26 21:23:09
% DurationCPUTime: 0.23s
% Computational Cost: add. (210->50), mult. (407->73), div. (0->0), fcn. (317->8), ass. (0->37)
t209 = sin(qJ(2));
t211 = cos(qJ(2));
t230 = t211 * qJD(5);
t237 = qJ(3) + cos(pkin(9)) * pkin(5) + pkin(4);
t227 = pkin(2) + pkin(3) + r_i_i_C(3) + pkin(8) + qJ(5);
t241 = t227 * t209;
t246 = (-t237 * t211 + t241) * qJD(2) + (pkin(5) * sin(pkin(9)) - pkin(7) + qJ(4)) * qJD(1) - t209 * qJD(3) - t230;
t206 = pkin(9) + qJ(6);
t204 = sin(t206);
t205 = cos(t206);
t218 = r_i_i_C(1) * t205 - r_i_i_C(2) * t204 + t237;
t244 = -t218 * t211 + t241;
t212 = cos(qJ(1));
t226 = qJD(6) * t209 + qJD(1);
t243 = t212 * t226;
t225 = qJD(1) * t209 + qJD(6);
t210 = sin(qJ(1));
t233 = qJD(2) * t211;
t229 = t210 * t233;
t239 = t225 * t212 + t229;
t236 = qJD(1) * t210;
t235 = qJD(1) * t212;
t234 = qJD(2) * t209;
t232 = qJD(2) * t212;
t231 = qJD(6) * t211;
t228 = t211 * t232;
t222 = t227 * t211;
t220 = t226 * t210;
t217 = qJD(3) + (-r_i_i_C(1) * t204 - r_i_i_C(2) * t205) * qJD(6);
t216 = t225 * t210 - t228;
t215 = -qJD(4) + (-t237 * t209 - pkin(1) - t222) * qJD(1);
t213 = -qJD(5) * t209 + t217 * t211 + (-t218 * t209 - t222) * qJD(2);
t202 = t204 * t220 - t239 * t205;
t201 = t239 * t204 + t205 * t220;
t200 = t204 * t243 + t216 * t205;
t199 = t216 * t204 - t205 * t243;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t246 * t210 + t215 * t212, t213 * t212 + t244 * t236, -t209 * t236 + t228, -t235, -t209 * t232 - t211 * t236, r_i_i_C(1) * t199 + r_i_i_C(2) * t200; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t215 * t210 - t246 * t212, t213 * t210 - t235 * t244, t209 * t235 + t229, -t236, -t210 * t234 + t211 * t235, -r_i_i_C(1) * t201 + r_i_i_C(2) * t202; 0, -qJD(2) * t244 + t217 * t209 + t230, t234, 0, t233 (-t204 * t231 - t205 * t234) * r_i_i_C(2) + (-t204 * t234 + t205 * t231) * r_i_i_C(1);];
JaD_transl  = t1;
