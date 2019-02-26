% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,d6,theta3]';
%
% Output:
% JaD_transl [3x6]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S6RRPPRR8_jacobiaD_transl_4_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_transl_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRR8_jacobiaD_transl_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S6RRPPRR8_jacobiaD_transl_4_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRR8_jacobiaD_transl_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:32:35
% EndTime: 2019-02-26 21:32:35
% DurationCPUTime: 0.21s
% Computational Cost: add. (88->38), mult. (286->64), div. (0->0), fcn. (229->6), ass. (0->32)
t179 = sin(qJ(2));
t177 = sin(pkin(10));
t178 = cos(pkin(10));
t203 = r_i_i_C(3) + qJ(4);
t206 = pkin(3) + r_i_i_C(1);
t207 = t177 * t203 + t178 * t206 + pkin(2);
t181 = cos(qJ(2));
t204 = r_i_i_C(2) + qJ(3);
t208 = t204 * t181;
t211 = t179 * t207 - t208;
t194 = qJD(4) * t177;
t187 = t179 * qJD(3) + t181 * t194;
t210 = (-pkin(2) * t179 + t208) * qJD(2) + t187;
t182 = cos(qJ(1));
t202 = t177 * t182;
t180 = sin(qJ(1));
t201 = t180 * t181;
t200 = t182 * t178;
t199 = qJD(1) * t180;
t198 = qJD(1) * t182;
t197 = qJD(2) * t179;
t196 = qJD(2) * t181;
t195 = qJD(2) * t182;
t193 = t178 * qJD(4);
t192 = t180 * t197;
t191 = t179 * t195;
t190 = t204 * t179;
t186 = -pkin(2) * t181 - pkin(1) - t190;
t183 = -t179 * t194 + qJD(3) * t181 + (-t181 * t207 - t190) * qJD(2);
t175 = -t177 * t192 + (-t178 * t180 + t181 * t202) * qJD(1);
t173 = t177 * t191 + (t177 * t201 + t200) * qJD(1);
t1 = [-t182 * t193 + t206 * (t178 * t192 + (-t177 * t180 - t181 * t200) * qJD(1)) - t203 * t175 - t210 * t180 + (-t180 * pkin(7) + t186 * t182) * qJD(1), t182 * t183 + t199 * t211, -t179 * t199 + t181 * t195, -t173, 0, 0; -t180 * t193 + t206 * (-t178 * t191 + (-t178 * t201 + t202) * qJD(1)) - t203 * t173 + t210 * t182 + (t182 * pkin(7) + t186 * t180) * qJD(1), t180 * t183 - t198 * t211, t179 * t198 + t180 * t196, t175, 0, 0; 0, -qJD(2) * t211 + t187, t197, t177 * t196, 0, 0;];
JaD_transl  = t1;
