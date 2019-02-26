% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPRPP5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:37
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRPP5_jacobiaD_transl_6_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_jacobiaD_transl_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRPP5_jacobiaD_transl_6_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_jacobiaD_transl_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:37:08
% EndTime: 2019-02-26 21:37:08
% DurationCPUTime: 0.32s
% Computational Cost: add. (229->55), mult. (694->85), div. (0->0), fcn. (577->6), ass. (0->39)
t205 = sin(qJ(2));
t208 = cos(qJ(2));
t207 = cos(qJ(4));
t222 = pkin(2) + pkin(8) - r_i_i_C(3) - qJ(6);
t240 = t222 * qJD(2) + qJD(5) * t207 - qJD(3);
t245 = t240 * t205 - (pkin(3) + pkin(7)) * qJD(1) - (qJ(3) * qJD(2) - qJD(6)) * t208;
t204 = sin(qJ(4));
t226 = r_i_i_C(1) + pkin(5) + pkin(4);
t237 = r_i_i_C(2) + qJ(5);
t214 = t226 * t204 - t237 * t207 + qJ(3);
t244 = -t222 * t205 + t214 * t208;
t243 = t214 * qJD(2) - qJD(6);
t206 = sin(qJ(1));
t229 = qJD(2) * t208;
t223 = t206 * t229;
t209 = cos(qJ(1));
t231 = qJD(1) * t209;
t241 = t205 * t231 + t223;
t236 = t204 * t209;
t235 = t206 * t204;
t234 = t206 * t207;
t233 = t207 * t209;
t232 = qJD(1) * t206;
t230 = qJD(2) * t205;
t228 = qJD(2) * t209;
t227 = qJD(4) * t208;
t224 = t208 * t228;
t220 = qJD(4) * t205 + qJD(1);
t219 = qJD(1) * t205 + qJD(4);
t217 = t220 * t207;
t216 = t205 * t236 + t234;
t212 = t204 * qJD(5) + (-qJ(3) * t205 - t222 * t208 - pkin(1)) * qJD(1);
t211 = (t237 * t204 + t226 * t207) * qJD(4) - t240;
t210 = -t243 * t205 + t211 * t208;
t199 = t209 * t217 + (-t219 * t206 + t224) * t204;
t198 = -t207 * t224 + t216 * qJD(4) + (t205 * t234 + t236) * qJD(1);
t197 = t206 * t217 + (t219 * t209 + t223) * t204;
t196 = -qJD(4) * t233 - t241 * t207 + t220 * t235;
t1 = [-t237 * t196 - t226 * t197 + t245 * t206 + t212 * t209, t210 * t209 - t244 * t232, -t205 * t232 + t224, t216 * qJD(5) - t226 * t198 + t237 * t199, t198, t205 * t228 + t208 * t232; t237 * t198 + t226 * t199 + t212 * t206 - t245 * t209, t210 * t206 + t244 * t231, t241 -(-t205 * t235 + t233) * qJD(5) + t237 * t197 - t226 * t196, t196, t206 * t230 - t208 * t231; 0, t211 * t205 + t243 * t208, t230 (t226 * t230 - t237 * t227) * t207 + (t237 * t230 + (t226 * qJD(4) - qJD(5)) * t208) * t204, -t204 * t227 - t207 * t230, -t229;];
JaD_transl  = t1;
