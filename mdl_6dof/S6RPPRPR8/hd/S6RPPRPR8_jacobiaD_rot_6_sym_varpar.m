% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRPR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:29
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRPR8_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:29:44
% EndTime: 2019-02-26 20:29:45
% DurationCPUTime: 0.67s
% Computational Cost: add. (1892->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
t137 = cos(qJ(1));
t133 = t137 ^ 2;
t131 = pkin(9) + qJ(4);
t129 = sin(t131);
t125 = t129 ^ 2;
t130 = cos(t131);
t127 = 0.1e1 / t130 ^ 2;
t182 = t125 * t127;
t121 = t133 * t182 + 0.1e1;
t124 = t129 * t125;
t126 = 0.1e1 / t130;
t180 = t126 * t129;
t145 = qJD(4) * (t124 * t126 * t127 + t180);
t135 = sin(qJ(1));
t172 = qJD(1) * t137;
t149 = t125 * t135 * t172;
t189 = (-t127 * t149 + t133 * t145) / t121 ^ 2;
t200 = -0.2e1 * t189;
t158 = 0.1e1 + t182;
t199 = t137 * t158;
t152 = qJD(1) * t130 + qJD(6);
t169 = qJD(4) * t137;
t198 = t129 * t169 + t152 * t135;
t174 = t137 * t129;
t120 = atan2(t174, t130);
t116 = sin(t120);
t117 = cos(t120);
t102 = t116 * t174 + t117 * t130;
t99 = 0.1e1 / t102;
t136 = cos(qJ(6));
t175 = t136 * t137;
t134 = sin(qJ(6));
t177 = t135 * t134;
t115 = -t130 * t177 + t175;
t109 = 0.1e1 / t115;
t100 = 0.1e1 / t102 ^ 2;
t110 = 0.1e1 / t115 ^ 2;
t118 = 0.1e1 / t121;
t197 = t118 - 0.1e1;
t132 = t135 ^ 2;
t171 = qJD(4) * t130;
t181 = t125 * t132;
t161 = t127 * t169;
t173 = qJD(1) * t135;
t162 = t129 * t173;
t93 = ((t130 * t169 - t162) * t126 + t125 * t161) * t118;
t154 = -t93 + t169;
t155 = t137 * t93 - qJD(4);
t184 = t117 * t129;
t88 = t155 * t184 + (t154 * t130 - t162) * t116;
t194 = t99 * t100 * t88;
t96 = t100 * t181 + 0.1e1;
t196 = (-t181 * t194 + (t129 * t132 * t171 + t149) * t100) / t96 ^ 2;
t94 = 0.1e1 / t96;
t195 = t100 * t94;
t176 = t135 * t136;
t178 = t134 * t137;
t114 = t130 * t176 + t178;
t108 = t114 ^ 2;
t107 = t108 * t110 + 0.1e1;
t187 = t110 * t114;
t153 = qJD(6) * t130 + qJD(1);
t147 = t153 * t136;
t170 = qJD(4) * t135;
t98 = -t135 * t147 + (t129 * t170 - t152 * t137) * t134;
t192 = t109 * t110 * t98;
t163 = t130 * t175;
t97 = -qJD(1) * t163 - qJD(6) * t175 + (qJD(4) * t129 * t136 + t153 * t134) * t135;
t193 = 0.1e1 / t107 ^ 2 * (-t108 * t192 - t97 * t187);
t188 = t109 * t136;
t186 = t114 * t134;
t185 = t116 * t130;
t183 = t125 * t126;
t179 = t129 * t135;
t168 = 0.2e1 * t196;
t167 = 0.2e1 * t194;
t166 = 0.2e1 * t193;
t165 = t94 * t171;
t164 = t137 * t183;
t159 = -0.2e1 * t99 * t196;
t157 = 0.2e1 * t114 * t192;
t156 = 0.2e1 * t129 * t189;
t151 = t118 * t164;
t150 = t197 * t129 * t116;
t148 = t135 * t158;
t146 = t110 * t186 + t188;
t113 = -t130 * t178 - t176;
t112 = t163 - t177;
t105 = 0.1e1 / t107;
t104 = t118 * t199;
t92 = (-t117 * t151 + t150) * t135;
t90 = t137 * t185 - t184 + (t117 * t174 - t185) * t104;
t89 = t199 * t200 + (-qJD(1) * t148 + 0.2e1 * t137 * t145) * t118;
t1 = [t135 * t126 * t156 + (-qJD(4) * t148 - t172 * t180) * t118, 0, 0, t89, 0, 0; (t99 * t165 + (t159 + (-qJD(1) * t92 - t88) * t195) * t129) * t137 + ((-t92 * t165 + (t92 * t168 + ((-t93 * t151 - t197 * t171 + t156) * t116 + (t164 * t200 + t129 * t93 + (t124 * t161 - (t93 - 0.2e1 * t169) * t129) * t118) * t117) * t94 * t135) * t129) * t100 + (t92 * t167 + (-t99 + ((-t132 + t133) * t118 * t117 * t183 - t137 * t150) * t100) * qJD(1)) * t129 * t94) * t135, 0, 0 (t99 * t94 * t172 + (t159 + (-qJD(4) * t90 - t88) * t195) * t135) * t130 + (((-qJD(4) * t99 + t90 * t167) * t135 + (-t90 * t172 + (-(-t104 * t173 + t137 * t89) * t117 - ((-t104 * t137 + 0.1e1) * t93 + (t104 - t137) * qJD(4)) * t116) * t179) * t100) * t94 + (t90 * t168 - ((-t89 - t173) * t116 + (t154 * t104 + t155) * t117) * t94 * t130) * t100 * t135) * t129, 0, 0; (-t109 * t112 + t113 * t187) * t166 + (t113 * t157 - t153 * t109 * t178 - t198 * t188 + (t137 * t114 * t147 - t112 * t98 + t113 * t97 - t198 * t186) * t110) * t105, 0, 0, t146 * t166 * t179 + (-t146 * t130 * t170 + (-t146 * t172 + ((qJD(6) * t109 + t157) * t134 + (t134 * t97 + (-qJD(6) * t114 + t98) * t136) * t110) * t135) * t129) * t105, 0, -0.2e1 * t193 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t192 - t110 * t193) * t114) * t114;];
JaD_rot  = t1;
