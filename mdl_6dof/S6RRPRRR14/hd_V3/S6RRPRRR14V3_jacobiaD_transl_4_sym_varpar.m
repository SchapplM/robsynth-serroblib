% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRRR14V3
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
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_jacobiaD_transl_4_sym_varpar: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:12:05
% EndTime: 2019-04-12 15:12:05
% DurationCPUTime: 0.12s
% Computational Cost: add. (77->35), mult. (254->70), div. (0->0), fcn. (208->6), ass. (0->31)
t197 = sin(qJ(2));
t196 = sin(qJ(4));
t199 = cos(qJ(4));
t207 = r_i_i_C(1) * t199 - r_i_i_C(2) * t196;
t200 = cos(qJ(2));
t222 = r_i_i_C(3) + qJ(3);
t223 = t222 * t200;
t203 = -t207 * t197 + t223;
t201 = cos(qJ(1));
t221 = t199 * t201;
t198 = sin(qJ(1));
t220 = qJD(1) * t198;
t219 = qJD(1) * t201;
t218 = qJD(2) * t197;
t217 = qJD(2) * t198;
t216 = qJD(2) * t200;
t215 = qJD(2) * t201;
t214 = qJD(4) * t197;
t213 = qJD(4) * t200;
t210 = qJD(1) * t222;
t209 = -qJD(1) + t213;
t208 = qJD(1) * t200 - qJD(4);
t206 = r_i_i_C(1) * t196 + r_i_i_C(2) * t199;
t205 = t209 * t196;
t204 = t197 * t215 + t208 * t198;
t202 = qJD(3) * t200 + t206 * t214 + (-t222 * t197 - t207 * t200) * qJD(2);
t195 = -t208 * t221 + (t199 * t218 + t205) * t198;
t194 = t209 * t199 * t198 + (-t197 * t217 + t208 * t201) * t196;
t193 = t204 * t199 + t201 * t205;
t192 = t204 * t196 - t209 * t221;
t1 = [t195 * r_i_i_C(1) + t194 * r_i_i_C(2) - t223 * t217 + (-qJD(3) * t198 - t201 * t210) * t197, t202 * t201 - t203 * t220, -t197 * t220 + t200 * t215, t192 * r_i_i_C(1) + t193 * r_i_i_C(2), 0, 0; -t193 * r_i_i_C(1) + t192 * r_i_i_C(2) + t223 * t215 + (qJD(3) * t201 - t198 * t210) * t197, t202 * t198 + t203 * t219, t197 * t219 + t198 * t216, -t194 * r_i_i_C(1) + t195 * r_i_i_C(2), 0, 0; 0, t203 * qJD(2) + t197 * qJD(3) - t206 * t213, t218 (t196 * t214 - t199 * t216) * r_i_i_C(2) + (-t196 * t216 - t199 * t214) * r_i_i_C(1), 0, 0;];
JaD_transl  = t1;
