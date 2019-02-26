% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRPPPR3
% Use Code from Maple symbolic Code Generation
%
% analytische Jacobi-Matrix: Differentieller Zusammenhang zwischen
% Endeffektorposition und verallgemeinerten Koordinaten.
% Zeitableitung der Winkeldarstellung des Endeffektors in Basis-Koordinaten
%
% Winkeldarstellung: Euler-XYZ-Winkel, rotx(alpha)*roty(beta)*rotz(gamma)
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:23
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPPR3_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR3_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR3_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:23:08
% EndTime: 2019-02-26 21:23:09
% DurationCPUTime: 0.61s
% Computational Cost: add. (1161->93), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->93)
t135 = sin(qJ(2));
t136 = sin(qJ(1));
t137 = cos(qJ(2));
t181 = t136 * t137;
t120 = atan2(-t181, t135);
t117 = cos(t120);
t116 = sin(t120);
t167 = t116 * t181;
t105 = t117 * t135 - t167;
t102 = 0.1e1 / t105;
t127 = pkin(9) + qJ(6);
t126 = cos(t127);
t138 = cos(qJ(1));
t184 = t135 * t138;
t165 = t126 * t184;
t125 = sin(t127);
t183 = t136 * t125;
t115 = t165 - t183;
t109 = 0.1e1 / t115;
t128 = 0.1e1 / t135;
t103 = 0.1e1 / t105 ^ 2;
t110 = 0.1e1 / t115 ^ 2;
t129 = 0.1e1 / t135 ^ 2;
t131 = t136 ^ 2;
t133 = t137 ^ 2;
t186 = t129 * t133;
t123 = t131 * t186 + 0.1e1;
t121 = 0.1e1 / t123;
t198 = t121 - 0.1e1;
t182 = t136 * t126;
t114 = t125 * t184 + t182;
t108 = t114 ^ 2;
t101 = t108 * t110 + 0.1e1;
t192 = t110 * t114;
t153 = qJD(1) * t135 + qJD(6);
t154 = qJD(6) * t135 + qJD(1);
t173 = qJD(2) * t138;
t162 = t137 * t173;
t188 = t125 * t138;
t94 = -t154 * t188 + (-t153 * t136 + t162) * t126;
t194 = t109 * t110 * t94;
t177 = qJD(1) * t138;
t93 = -qJD(6) * t165 - t125 * t162 - t126 * t177 + t153 * t183;
t197 = (-t108 * t194 - t93 * t192) / t101 ^ 2;
t134 = t138 ^ 2;
t185 = t133 * t134;
t100 = t103 * t185 + 0.1e1;
t96 = 0.1e1 / t100;
t196 = t103 * t96;
t175 = qJD(2) * t136;
t163 = t129 * t175;
t164 = t137 * t177;
t95 = ((t135 * t175 - t164) * t128 + t133 * t163) * t121;
t155 = -t95 + t175;
t156 = -t136 * t95 + qJD(2);
t189 = t117 * t137;
t89 = t156 * t189 + (t155 * t135 - t164) * t116;
t195 = t102 * t103 * t89;
t193 = t109 * t125;
t191 = t114 * t126;
t190 = t116 * t135;
t187 = t128 * t133;
t179 = qJD(1) * t136;
t178 = qJD(1) * t137;
t176 = qJD(2) * t135;
t174 = qJD(2) * t137;
t151 = t133 * t136 * t177;
t172 = 0.2e1 * (-t185 * t195 + (-t134 * t135 * t174 - t151) * t103) / t100 ^ 2;
t171 = 0.2e1 * t197;
t170 = 0.2e1 * t195;
t132 = t137 * t133;
t148 = qJD(2) * (-t129 * t132 - t137) * t128;
t169 = 0.2e1 * (t129 * t151 + t131 * t148) / t123 ^ 2;
t168 = t96 * t176;
t166 = t121 * t187;
t161 = 0.1e1 + t186;
t160 = t102 * t172;
t159 = 0.2e1 * t114 * t194;
t158 = t136 * t169;
t157 = t137 * t169;
t152 = t136 * t166;
t150 = t161 * t138;
t149 = t110 * t191 - t193;
t98 = 0.1e1 / t101;
t147 = t149 * t98;
t146 = -t136 * t174 - t153 * t138;
t113 = -t135 * t182 - t188;
t112 = t126 * t138 - t135 * t183;
t107 = t161 * t136 * t121;
t92 = (t198 * t137 * t116 + t117 * t152) * t138;
t91 = t136 * t190 + t189 + (-t117 * t181 - t190) * t107;
t90 = -t161 * t158 + (qJD(1) * t150 + 0.2e1 * t136 * t148) * t121;
t1 = [t128 * t138 * t157 + (t128 * t136 * t178 + qJD(2) * t150) * t121, t90, 0, 0, 0, 0; (t102 * t168 + (t160 + (qJD(1) * t92 + t89) * t196) * t137) * t136 + (t92 * t137 * t96 * t170 + (t92 * t168 + (t92 * t172 + ((t95 * t152 + t198 * t176 + t157) * t116 + (t158 * t187 + t137 * t95 + (t132 * t163 - (t95 - 0.2e1 * t175) * t137) * t121) * t117) * t96 * t138) * t137) * t103 + (-t102 + ((t131 - t134) * t117 * t166 + t198 * t167) * t103) * t96 * t178) * t138 (t102 * t96 * t179 + (t160 + (qJD(2) * t91 + t89) * t196) * t138) * t135 + (t91 * t138 * t103 * t172 + ((-qJD(2) * t102 + t91 * t170) * t138 + (t91 * t179 + (-(-t107 * t177 - t136 * t90) * t117 - ((t107 * t136 - 0.1e1) * t95 + (-t107 + t136) * qJD(2)) * t116) * t137 * t138) * t103) * t96 - ((-t90 + t177) * t116 + (t155 * t107 - t156) * t117) * t184 * t196) * t137, 0, 0, 0, 0; (-t109 * t112 + t113 * t192) * t171 + (t113 * t159 - t154 * t109 * t182 + t146 * t193 + (-t154 * t114 * t183 - t112 * t94 + t113 * t93 - t146 * t191) * t110) * t98, t135 * t147 * t173 + (t147 * t179 + (t149 * t171 + ((qJD(6) * t109 + t159) * t126 + (t126 * t93 + (qJD(6) * t114 - t94) * t125) * t110) * t98) * t138) * t137, 0, 0, 0, -0.2e1 * t197 + 0.2e1 * (-t110 * t93 * t98 + (-t110 * t197 - t98 * t194) * t114) * t114;];
JaD_rot  = t1;
