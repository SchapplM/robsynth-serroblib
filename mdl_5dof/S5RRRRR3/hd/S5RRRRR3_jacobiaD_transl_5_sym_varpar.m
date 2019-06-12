% Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% r_i_i_C [3x1]
%   Ortsvektor vom KörperKS-Ursprung zum gesuchten Punkt
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
%
% Output:
% JaD_transl [3x5]
%   Zeitableitung der translatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-31 10:31
% Revision: 36f6366a01c4a552c0708fcd8ed3e0fb9da693e2 (2019-05-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_transl = S5RRRRR3_jacobiaD_transl_5_sym_varpar(qJ, qJD, r_i_i_C, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(3,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: qJD has to be [5x1] (double)');
assert(isa(r_i_i_C,'double') && isreal(r_i_i_C) && all(size(r_i_i_C) == [3 1]), ...
	'S5RRRRR3_jacobiaD_transl_5_sym_varpar: Position vector r_i_i_C has to be [3x1] double');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_jacobiaD_transl_5_sym_varpar: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From jacobiaD_transl_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-31 10:31:14
% EndTime: 2019-05-31 10:31:14
% DurationCPUTime: 0.28s
% Computational Cost: add. (504->72), mult. (556->111), div. (0->0), fcn. (429->10), ass. (0->70)
t133 = qJ(2) + qJ(3);
t127 = sin(t133);
t137 = cos(qJ(4));
t125 = t137 * pkin(3) + pkin(2);
t132 = qJ(4) + qJ(5);
t128 = cos(t132);
t196 = r_i_i_C(1) * t128 + t125;
t200 = t127 * t196;
t126 = sin(t132);
t129 = cos(t133);
t131 = qJD(2) + qJD(3);
t180 = t129 * t131;
t130 = qJD(4) + qJD(5);
t181 = t128 * t130;
t199 = t126 * t180 + t127 * t181;
t191 = pkin(5) + r_i_i_C(3);
t198 = t191 * t129;
t168 = t191 * t131;
t135 = sin(qJ(2));
t186 = pkin(1) * qJD(2);
t170 = t135 * t186;
t134 = sin(qJ(4));
t185 = pkin(3) * qJD(4);
t173 = t134 * t185;
t183 = t127 * t131;
t189 = pkin(3) * t134;
t197 = qJD(1) * t189 - t125 * t183 + (t168 - t173) * t129 - t170;
t184 = t126 * t127;
t166 = t130 * t184;
t194 = r_i_i_C(1) * t166 + t199 * r_i_i_C(2) + t127 * t173;
t139 = cos(qJ(1));
t155 = t129 * t130 - qJD(1);
t193 = t139 * t155;
t177 = qJD(1) * t129;
t154 = -t130 + t177;
t136 = sin(qJ(1));
t163 = t136 * t183;
t192 = t154 * t139 - t163;
t190 = pkin(1) * t135;
t187 = r_i_i_C(2) * t128;
t182 = t127 * t139;
t162 = t131 * t182;
t144 = t154 * t136 + t162;
t104 = t144 * t126 - t128 * t193;
t105 = t126 * t193 + t144 * t128;
t179 = t104 * r_i_i_C(1) + t105 * r_i_i_C(2);
t148 = t155 * t136;
t106 = t192 * t126 + t128 * t148;
t107 = t126 * t148 - t192 * t128;
t178 = -t106 * r_i_i_C(1) + t107 * r_i_i_C(2);
t176 = qJD(1) * t136;
t175 = qJD(1) * t139;
t174 = r_i_i_C(2) * t184;
t172 = t137 * t185;
t171 = r_i_i_C(2) * qJD(1) * t126;
t169 = t191 * t127;
t167 = t191 * t136;
t152 = -qJD(4) + t177;
t151 = -r_i_i_C(1) * t126 - t187;
t150 = t196 * t131;
t149 = t196 * t139;
t147 = (-qJD(4) * t129 + qJD(1)) * t137;
t146 = t194 * t139 + t176 * t200;
t145 = t194 * t136 + t171 * t182 + t175 * t198;
t138 = cos(qJ(2));
t142 = t172 + (-pkin(1) * t138 - t125 * t129 - t169) * qJD(1);
t141 = -t138 * t186 + (-t129 * t196 - t169) * t131;
t140 = t131 * t174 - t127 * t150 + (t151 * t130 - t173) * t129 + t191 * t180;
t117 = r_i_i_C(2) * t166;
t1 = [t107 * r_i_i_C(1) + t106 * r_i_i_C(2) - t197 * t136 + t142 * t139, (-t174 + t190 - t198) * t176 + t141 * t139 + t146, (-t136 * t171 - t139 * t168) * t127 + (-qJD(1) * t167 - t131 * t149) * t129 + t146, (t139 * t147 + (t152 * t136 + t162) * t134) * pkin(3) + t179, t179; -t105 * r_i_i_C(1) + t104 * r_i_i_C(2) + t142 * t136 + t197 * t139, (-t190 - t200) * t175 + t141 * t136 + t145, -t136 * t129 * t150 + (-qJD(1) * t149 - t131 * t167) * t127 + t145, (t136 * t147 + (-t152 * t139 + t163) * t134) * pkin(3) + t178, t178; 0, t140 - t170, t140, t117 + (-r_i_i_C(1) * t181 - t172) * t127 + (t151 - t189) * t180, -t199 * r_i_i_C(1) - t180 * t187 + t117;];
JaD_transl  = t1;
