% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:26
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRP3_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP3_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRP3_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPPRP3_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:26:10
% EndTime: 2019-02-26 21:26:10
% DurationCPUTime: 0.27s
% Computational Cost: add. (169->52), mult. (474->74), div. (0->0), fcn. (365->6), ass. (0->39)
t240 = pkin(5) + r_i_i_C(1);
t205 = sin(qJ(5));
t207 = sin(qJ(1));
t210 = cos(qJ(1));
t206 = sin(qJ(2));
t226 = qJD(1) * t206 + qJD(5);
t218 = t226 * t210;
t208 = cos(qJ(5));
t227 = qJD(5) * t206 + qJD(1);
t220 = t227 * t208;
t209 = cos(qJ(2));
t233 = qJD(2) * t209;
t230 = t207 * t233;
t201 = (t218 + t230) * t205 + t207 * t220;
t231 = t209 * qJD(6);
t237 = t208 * pkin(5) + pkin(4) + qJ(3);
t238 = pkin(5) * qJD(5);
t228 = pkin(2) + pkin(3) + r_i_i_C(3) + qJ(6) + pkin(8);
t243 = t228 * t206;
t249 = qJD(2) * (-t237 * t209 + t243) + (pkin(5) * t205 - pkin(7) + qJ(4)) * qJD(1) + (t205 * t238 - qJD(3)) * t206 - t231;
t239 = r_i_i_C(2) * t205;
t216 = r_i_i_C(1) * t208 + t237 - t239;
t246 = -t216 * t209 + t243;
t236 = qJD(1) * t207;
t235 = qJD(1) * t210;
t234 = qJD(2) * t206;
t232 = qJD(2) * t210;
t229 = t209 * t232;
t222 = t228 * t209;
t219 = t227 * t205;
t217 = -r_i_i_C(2) * t208 - t240 * t205;
t215 = t217 * qJD(5) + qJD(3);
t214 = t226 * t207 - t229;
t212 = -t208 * t238 - qJD(4) + (-t237 * t206 - pkin(1) - t222) * qJD(1);
t199 = t214 * t205 - t210 * t220;
t211 = -qJD(6) * t206 + t215 * t209 + (-t216 * t206 - t222) * qJD(2);
t202 = -t208 * t218 + (-t208 * t233 + t219) * t207;
t200 = t214 * t208 + t210 * t219;
t1 = [t202 * r_i_i_C(1) + t201 * r_i_i_C(2) + t249 * t207 + t212 * t210, t211 * t210 + t236 * t246, -t206 * t236 + t229, -t235, t200 * r_i_i_C(2) + t240 * t199, -t206 * t232 - t209 * t236; -t200 * r_i_i_C(1) + t199 * r_i_i_C(2) + t212 * t207 - t249 * t210, t211 * t207 - t235 * t246, t206 * t235 + t230, -t236, t202 * r_i_i_C(2) - t240 * t201, -t207 * t234 + t209 * t235; 0, -qJD(2) * t246 + t215 * t206 + t231, t234, 0 (t240 * t208 - t239) * t209 * qJD(5) + t217 * t234, t233;];
JaD_transl  = t1;
