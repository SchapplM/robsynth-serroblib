% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta4]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR4_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR4_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR4_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR4_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:36
% EndTime: 2019-02-26 21:23:36
% DurationCPUTime: 0.18s
% Computational Cost: add. (107->41), mult. (334->63), div. (0->0), fcn. (266->6), ass. (0->30)
t176 = sin(qJ(2));
t178 = cos(qJ(2));
t189 = pkin(2) + r_i_i_C(2) + qJ(4);
t174 = sin(pkin(9));
t175 = cos(pkin(9));
t199 = r_i_i_C(3) + qJ(5);
t200 = pkin(4) + r_i_i_C(1);
t202 = t200 * t174 - t199 * t175 + qJ(3);
t181 = -t189 * t176 + t202 * t178;
t185 = -qJD(5) * t175 + qJD(3);
t206 = (t189 * qJD(2) - t185) * t176 - (pkin(3) + pkin(7)) * qJD(1) - (qJ(3) * qJD(2) + qJD(4)) * t178;
t177 = sin(qJ(1));
t198 = t175 * t177;
t179 = cos(qJ(1));
t197 = t175 * t179;
t196 = t177 * t174;
t195 = t179 * t174;
t194 = qJD(1) * t177;
t193 = qJD(1) * t179;
t192 = qJD(2) * t176;
t191 = qJD(2) * t178;
t190 = qJD(2) * t179;
t188 = t177 * t191;
t187 = t178 * t190;
t184 = t189 * t178;
t182 = t174 * qJD(5) + (-qJ(3) * t176 - pkin(1) - t184) * qJD(1);
t180 = -qJD(4) * t176 + t185 * t178 + (-t176 * t202 - t184) * qJD(2);
t172 = -t175 * t187 + (t176 * t198 + t195) * qJD(1);
t170 = -t175 * t188 + (-t176 * t197 + t196) * qJD(1);
t1 = [t200 * (-t174 * t188 + (-t176 * t195 - t198) * qJD(1)) - t199 * t170 + t182 * t179 + t206 * t177, t180 * t179 - t181 * t194, -t176 * t194 + t187, -t176 * t190 - t178 * t194, t172, 0; t200 * (t174 * t187 + (-t176 * t196 + t197) * qJD(1)) + t199 * t172 + t182 * t177 - t206 * t179, t180 * t177 + t181 * t193, t176 * t193 + t188, -t177 * t192 + t178 * t193, t170, 0; 0, t181 * qJD(2) + t178 * qJD(4) + t185 * t176, t192, t191, -t175 * t192, 0;];
JaD_transl  = t1;
