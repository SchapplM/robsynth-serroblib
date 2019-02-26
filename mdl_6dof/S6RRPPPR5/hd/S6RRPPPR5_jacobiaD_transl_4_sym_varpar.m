% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
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

function JaD_transl = S6RRPPPR5_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPPR5_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_jacobiaD_transl_4_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:13
% EndTime: 2019-02-26 21:24:13
% DurationCPUTime: 0.16s
% Computational Cost: add. (88->38), mult. (286->64), div. (0->0), fcn. (229->6), ass. (0->32)
t184 = sin(qJ(2));
t182 = sin(pkin(9));
t183 = cos(pkin(9));
t208 = r_i_i_C(3) + qJ(4);
t211 = pkin(3) - r_i_i_C(2);
t212 = t208 * t182 + t211 * t183 + pkin(2);
t186 = cos(qJ(2));
t209 = r_i_i_C(1) + qJ(3);
t213 = t209 * t186;
t216 = t212 * t184 - t213;
t199 = qJD(4) * t182;
t192 = t184 * qJD(3) + t186 * t199;
t215 = (-pkin(2) * t184 + t213) * qJD(2) + t192;
t187 = cos(qJ(1));
t207 = t182 * t187;
t185 = sin(qJ(1));
t206 = t185 * t186;
t205 = t187 * t183;
t204 = qJD(1) * t185;
t203 = qJD(1) * t187;
t202 = qJD(2) * t184;
t201 = qJD(2) * t186;
t200 = qJD(2) * t187;
t198 = t183 * qJD(4);
t197 = t185 * t202;
t196 = t184 * t200;
t195 = t209 * t184;
t191 = -pkin(2) * t186 - pkin(1) - t195;
t188 = -t184 * t199 + qJD(3) * t186 + (-t186 * t212 - t195) * qJD(2);
t180 = -t182 * t197 + (-t185 * t183 + t186 * t207) * qJD(1);
t178 = t182 * t196 + (t182 * t206 + t205) * qJD(1);
t1 = [-t187 * t198 - t211 * (-t183 * t197 + (t182 * t185 + t186 * t205) * qJD(1)) - t208 * t180 - t215 * t185 + (-t185 * pkin(7) + t191 * t187) * qJD(1), t188 * t187 + t216 * t204, -t184 * t204 + t186 * t200, -t178, 0, 0; -t185 * t198 - t211 * (t183 * t196 + (t183 * t206 - t207) * qJD(1)) - t208 * t178 + t215 * t187 + (t187 * pkin(7) + t191 * t185) * qJD(1), t188 * t185 - t203 * t216, t184 * t203 + t185 * t201, t180, 0, 0; 0, -qJD(2) * t216 + t192, t202, t182 * t201, 0, 0;];
JaD_transl  = t1;
