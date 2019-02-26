% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPRPR8
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d6,theta3]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:41
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPRPR8_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPR8_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRPR8_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:41:47
% EndTime: 2019-02-26 21:41:48
% DurationCPUTime: 0.69s
% Computational Cost: add. (1350->93), mult. (2519->209), div. (480->12), fcn. (2968->9), ass. (0->94)
t133 = sin(qJ(1));
t127 = t133 ^ 2;
t132 = sin(qJ(2));
t126 = t132 ^ 2;
t134 = cos(qJ(2));
t129 = 0.1e1 / t134 ^ 2;
t182 = t126 * t129;
t119 = t127 * t182 + 0.1e1;
t125 = t132 * t126;
t128 = 0.1e1 / t134;
t179 = t128 * t132;
t144 = qJD(2) * (t125 * t128 * t129 + t179);
t135 = cos(qJ(1));
t171 = qJD(1) * t135;
t180 = t126 * t133;
t148 = t171 * t180;
t188 = (t127 * t144 + t129 * t148) / t119 ^ 2;
t198 = -0.2e1 * t188;
t124 = pkin(10) + qJ(4);
t123 = cos(t124);
t173 = t134 * t135;
t122 = sin(t124);
t177 = t133 * t122;
t112 = t123 * t173 + t177;
t107 = 0.1e1 / t112 ^ 2;
t176 = t133 * t123;
t111 = t122 * t173 - t176;
t185 = t111 * t123;
t106 = 0.1e1 / t112;
t187 = t106 * t122;
t146 = t107 * t185 - t187;
t105 = t111 ^ 2;
t98 = t105 * t107 + 0.1e1;
t96 = 0.1e1 / t98;
t197 = t146 * t96;
t155 = 0.1e1 + t182;
t196 = t133 * t155;
t175 = t133 * t132;
t116 = atan2(-t175, -t134);
t114 = cos(t116);
t113 = sin(t116);
t161 = t113 * t175;
t102 = -t114 * t134 - t161;
t99 = 0.1e1 / t102;
t100 = 0.1e1 / t102 ^ 2;
t117 = 0.1e1 / t119;
t195 = t117 - 0.1e1;
t131 = t135 ^ 2;
t169 = qJD(2) * t134;
t181 = t126 * t131;
t170 = qJD(2) * t133;
t157 = t129 * t170;
t158 = t132 * t171;
t92 = (-(-t133 * t169 - t158) * t128 + t126 * t157) * t117;
t152 = t92 - t170;
t153 = -t133 * t92 + qJD(2);
t184 = t114 * t132;
t86 = t153 * t184 + (t152 * t134 - t158) * t113;
t191 = t99 * t100 * t86;
t95 = t100 * t181 + 0.1e1;
t194 = (-t181 * t191 + (t131 * t132 * t169 - t148) * t100) / t95 ^ 2;
t186 = t107 * t111;
t150 = -qJD(1) * t134 + qJD(4);
t151 = qJD(4) * t134 - qJD(1);
t168 = qJD(2) * t135;
t156 = t132 * t168;
t183 = t122 * t135;
t91 = -t151 * t183 + (t150 * t133 - t156) * t123;
t190 = t106 * t107 * t91;
t174 = t133 * t134;
t145 = t122 * t174 + t123 * t135;
t90 = t145 * qJD(1) - t112 * qJD(4) + t122 * t156;
t193 = (-t105 * t190 - t90 * t186) / t98 ^ 2;
t93 = 0.1e1 / t95;
t192 = t100 * t93;
t178 = t132 * t135;
t172 = qJD(1) * t133;
t167 = 0.2e1 * t194;
t166 = -0.2e1 * t193;
t165 = 0.2e1 * t191;
t164 = t99 * t194;
t163 = t111 * t190;
t162 = t93 * t169;
t160 = t117 * t126 * t128;
t154 = t128 * t198;
t149 = t133 * t160;
t147 = t155 * t135;
t143 = t132 * t170 + t150 * t135;
t110 = -t123 * t174 + t183;
t104 = t117 * t196;
t89 = (t195 * t132 * t113 - t114 * t149) * t135;
t88 = -t113 * t174 + t184 + (t113 * t134 - t114 * t175) * t104;
t87 = t196 * t198 + (qJD(1) * t147 + 0.2e1 * t133 * t144) * t117;
t1 = [t154 * t178 + (qJD(2) * t147 - t172 * t179) * t117, t87, 0, 0, 0, 0; (-t99 * t162 + (0.2e1 * t164 + (qJD(1) * t89 + t86) * t192) * t132) * t133 + ((-t89 * t162 + (t89 * t167 + ((0.2e1 * t132 * t188 - t92 * t149 - t195 * t169) * t113 + (t154 * t180 + t132 * t92 + (t125 * t157 - (t92 - 0.2e1 * t170) * t132) * t117) * t114) * t93 * t135) * t132) * t100 + (t89 * t165 + (-t99 + ((-t127 + t131) * t114 * t160 + t195 * t161) * t100) * qJD(1)) * t132 * t93) * t135 (-t99 * t93 * t172 + (-0.2e1 * t164 + (-qJD(2) * t88 - t86) * t192) * t135) * t134 + (t88 * t135 * t100 * t167 + ((-qJD(2) * t99 + t88 * t165) * t135 + (t88 * t172 + (-(-t104 * t171 - t133 * t87) * t114 - ((t104 * t133 - 0.1e1) * t92 + (-t104 + t133) * qJD(2)) * t113) * t178) * t100) * t93 - ((t87 - t171) * t113 + (t152 * t104 + t153) * t114) * t173 * t192) * t132, 0, 0, 0, 0; 0.2e1 * (t106 * t145 + t110 * t186) * t193 + (0.2e1 * t110 * t163 - t151 * t106 * t176 + t143 * t187 + (-t151 * t111 * t177 + t110 * t90 - t143 * t185 + t145 * t91) * t107) * t96, t134 * t168 * t197 + (-t172 * t197 + (t146 * t166 + ((-qJD(4) * t106 - 0.2e1 * t163) * t123 + (-t123 * t90 + (-qJD(4) * t111 + t91) * t122) * t107) * t96) * t135) * t132, 0, t166 + 0.2e1 * (-t107 * t90 * t96 + (-t107 * t193 - t96 * t190) * t111) * t111, 0, 0;];
JaD_rot  = t1;
