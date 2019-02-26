% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3,theta5]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:34
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP1_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP1_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP1_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPP1_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:34:52
% EndTime: 2019-02-26 21:34:53
% DurationCPUTime: 0.16s
% Computational Cost: add. (167->41), mult. (296->69), div. (0->0), fcn. (230->8), ass. (0->34)
t212 = qJ(2) + pkin(9);
t210 = sin(t212);
t211 = cos(t212);
t243 = pkin(8) + r_i_i_C(3);
t225 = t243 * t211 - sin(qJ(2)) * pkin(2);
t246 = -pkin(3) * t210 + t225;
t214 = sin(qJ(4));
t217 = cos(qJ(4));
t227 = r_i_i_C(1) * t217 - r_i_i_C(2) * t214 + pkin(3);
t221 = -t227 * t210 + t225;
t244 = -t243 * t210 - cos(qJ(2)) * pkin(2);
t219 = cos(qJ(1));
t239 = t217 * t219;
t216 = sin(qJ(1));
t238 = qJD(1) * t216;
t237 = qJD(1) * t219;
t236 = qJD(2) * t216;
t235 = qJD(2) * t217;
t234 = qJD(2) * t219;
t233 = qJD(4) * t210;
t232 = qJD(4) * t211;
t230 = -qJD(1) + t232;
t229 = qJD(1) * t211 - qJD(4);
t228 = r_i_i_C(1) * t214 + r_i_i_C(2) * t217;
t226 = t230 * t214;
t223 = -pkin(3) * t211 - pkin(1) + t244;
t222 = t210 * t234 + t229 * t216;
t220 = t228 * t233 + (-t227 * t211 + t244) * qJD(2);
t213 = -qJ(3) - pkin(7);
t208 = -t229 * t239 + (t210 * t235 + t226) * t216;
t207 = t230 * t217 * t216 + (-t210 * t236 + t229 * t219) * t214;
t206 = t222 * t217 + t219 * t226;
t205 = t222 * t214 - t230 * t239;
t1 = [t208 * r_i_i_C(1) + t207 * r_i_i_C(2) + t219 * qJD(3) - t246 * t236 + (t213 * t216 + t223 * t219) * qJD(1), t220 * t219 - t221 * t238, t237, t205 * r_i_i_C(1) + t206 * r_i_i_C(2), 0, 0; -t206 * r_i_i_C(1) + t205 * r_i_i_C(2) + t216 * qJD(3) + t246 * t234 + (-t213 * t219 + t223 * t216) * qJD(1), t220 * t216 + t221 * t237, t238, -t207 * r_i_i_C(1) + t208 * r_i_i_C(2), 0, 0; 0, t221 * qJD(2) - t228 * t232, 0 (-t211 * t235 + t214 * t233) * r_i_i_C(2) + (-qJD(2) * t211 * t214 - t217 * t233) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
