% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRPPPR5
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPPR5_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_transl_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.18s
% Computational Cost: add. (119->42), mult. (383->68), div. (0->0), fcn. (312->6), ass. (0->31)
t180 = sin(qJ(2));
t182 = cos(qJ(2));
t195 = pkin(4) - r_i_i_C(2) + qJ(3);
t178 = sin(pkin(9));
t179 = cos(pkin(9));
t194 = -r_i_i_C(3) - qJ(5) - pkin(3);
t203 = r_i_i_C(1) + qJ(4);
t204 = t203 * t178 - t194 * t179 + pkin(2);
t209 = t204 * t180 - t195 * t182;
t188 = qJD(4) * t178 + qJD(5) * t179;
t208 = (t195 * qJD(2) + t188) * t182 - (pkin(2) * qJD(2) - qJD(3)) * t180;
t181 = sin(qJ(1));
t202 = t181 * t178;
t183 = cos(qJ(1));
t201 = t183 * t179;
t200 = qJD(1) * t181;
t199 = qJD(1) * t183;
t198 = qJD(2) * t180;
t197 = qJD(2) * t182;
t196 = qJD(2) * t183;
t193 = t181 * t198;
t192 = t180 * t196;
t190 = t195 * t180;
t189 = -t179 * qJD(4) + t178 * qJD(5);
t187 = -pkin(2) * t182 - pkin(1) - t190;
t184 = qJD(3) * t182 - t188 * t180 + (-t182 * t204 - t190) * qJD(2);
t175 = -t179 * t193 + (t182 * t201 + t202) * qJD(1);
t174 = -t178 * t193 + (t183 * t182 * t178 - t181 * t179) * qJD(1);
t173 = -t178 * t199 + (t182 * t200 + t192) * t179;
t172 = t178 * t192 + (t182 * t202 + t201) * qJD(1);
t1 = [t189 * t183 - t203 * t174 + t194 * t175 - t208 * t181 + (-t181 * pkin(7) + t187 * t183) * qJD(1), t184 * t183 + t209 * t200, -t180 * t200 + t182 * t196, -t172, -t173, 0; t189 * t181 - t203 * t172 + t194 * t173 + t208 * t183 + (t183 * pkin(7) + t187 * t181) * qJD(1), t184 * t181 - t199 * t209, t180 * t199 + t181 * t197, t174, t175, 0; 0, -qJD(2) * t209 + t180 * qJD(3) + t188 * t182, t198, t178 * t197, t179 * t197, 0;];
JaD_transl  = t1;
