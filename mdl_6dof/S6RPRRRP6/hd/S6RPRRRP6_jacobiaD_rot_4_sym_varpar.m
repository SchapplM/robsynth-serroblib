% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RPRRRP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:10
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRRP6_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP6_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP6_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:10:50
% EndTime: 2019-02-26 21:10:51
% DurationCPUTime: 0.72s
% Computational Cost: add. (2270->94), mult. (2519->210), div. (480->12), fcn. (2968->9), ass. (0->93)
t136 = sin(qJ(1));
t133 = t136 ^ 2;
t132 = pkin(10) + qJ(3);
t130 = sin(t132);
t126 = t130 ^ 2;
t131 = cos(t132);
t128 = 0.1e1 / t131 ^ 2;
t185 = t126 * t128;
t121 = t133 * t185 + 0.1e1;
t125 = t130 * t126;
t127 = 0.1e1 / t131;
t182 = t127 * t130;
t146 = qJD(3) * (t125 * t127 * t128 + t182);
t138 = cos(qJ(1));
t174 = qJD(1) * t138;
t183 = t126 * t136;
t151 = t174 * t183;
t191 = (t128 * t151 + t133 * t146) / t121 ^ 2;
t201 = -0.2e1 * t191;
t158 = 0.1e1 + t185;
t200 = t136 * t158;
t137 = cos(qJ(4));
t176 = t137 * t138;
t135 = sin(qJ(4));
t178 = t136 * t135;
t117 = t131 * t176 + t178;
t179 = t136 * t130;
t118 = atan2(-t179, -t131);
t113 = cos(t118);
t112 = sin(t118);
t164 = t112 * t179;
t102 = -t113 * t131 - t164;
t99 = 0.1e1 / t102;
t109 = 0.1e1 / t117;
t100 = 0.1e1 / t102 ^ 2;
t110 = 0.1e1 / t117 ^ 2;
t119 = 0.1e1 / t121;
t199 = t119 - 0.1e1;
t134 = t138 ^ 2;
t173 = qJD(3) * t131;
t184 = t126 * t134;
t172 = qJD(3) * t136;
t160 = t128 * t172;
t161 = t130 * t174;
t93 = (-(-t131 * t172 - t161) * t127 + t126 * t160) * t119;
t155 = t93 - t172;
t156 = -t136 * t93 + qJD(3);
t187 = t113 * t130;
t88 = t156 * t187 + (t131 * t155 - t161) * t112;
t196 = t99 * t100 * t88;
t96 = t100 * t184 + 0.1e1;
t198 = (-t184 * t196 + (t130 * t134 * t173 - t151) * t100) / t96 ^ 2;
t94 = 0.1e1 / t96;
t197 = t100 * t94;
t177 = t136 * t137;
t180 = t135 * t138;
t116 = t131 * t180 - t177;
t108 = t116 ^ 2;
t107 = t108 * t110 + 0.1e1;
t189 = t110 * t116;
t153 = -qJD(1) * t131 + qJD(4);
t154 = qJD(4) * t131 - qJD(1);
t171 = qJD(3) * t138;
t159 = t130 * t171;
t98 = -t154 * t180 + (t136 * t153 - t159) * t137;
t194 = t109 * t110 * t98;
t147 = t131 * t178 + t176;
t97 = t147 * qJD(1) - t117 * qJD(4) + t135 * t159;
t195 = 0.1e1 / t107 ^ 2 * (-t108 * t194 - t189 * t97);
t190 = t109 * t135;
t188 = t112 * t131;
t186 = t116 * t137;
t181 = t130 * t138;
t175 = qJD(1) * t136;
t170 = 0.2e1 * t198;
t169 = 0.2e1 * t196;
t168 = -0.2e1 * t195;
t167 = t99 * t198;
t166 = t116 * t194;
t165 = t94 * t173;
t163 = t119 * t126 * t127;
t157 = t127 * t201;
t152 = t136 * t163;
t150 = t158 * t138;
t149 = t153 * t138;
t148 = t110 * t186 - t190;
t115 = -t131 * t177 + t180;
t105 = 0.1e1 / t107;
t104 = t119 * t200;
t92 = (t112 * t130 * t199 - t113 * t152) * t138;
t90 = -t136 * t188 + t187 + (-t113 * t179 + t188) * t104;
t89 = t200 * t201 + (qJD(1) * t150 + 0.2e1 * t136 * t146) * t119;
t1 = [t157 * t181 + (qJD(3) * t150 - t175 * t182) * t119, 0, t89, 0, 0, 0; (-t99 * t165 + (0.2e1 * t167 + (qJD(1) * t92 + t88) * t197) * t130) * t136 + ((-t92 * t165 + (t92 * t170 + ((0.2e1 * t130 * t191 - t152 * t93 - t173 * t199) * t112 + (t157 * t183 + t130 * t93 + (t125 * t160 - (t93 - 0.2e1 * t172) * t130) * t119) * t113) * t94 * t138) * t130) * t100 + (t92 * t169 + (-t99 + ((-t133 + t134) * t113 * t163 + t199 * t164) * t100) * qJD(1)) * t130 * t94) * t138, 0 (-t99 * t94 * t175 + (-0.2e1 * t167 + (-qJD(3) * t90 - t88) * t197) * t138) * t131 + (((-qJD(3) * t99 + t169 * t90) * t138 + (t90 * t175 + (-(-t104 * t174 - t136 * t89) * t113 - ((t104 * t136 - 0.1e1) * t93 + (-t104 + t136) * qJD(3)) * t112) * t181) * t100) * t94 + (t90 * t170 - ((t89 - t174) * t112 + (t104 * t155 + t156) * t113) * t94 * t131) * t100 * t138) * t130, 0, 0, 0; 0.2e1 * (t109 * t147 + t115 * t189) * t195 + (0.2e1 * t115 * t166 - t154 * t109 * t177 + (t130 * t172 + t149) * t190 + (t147 * t98 + t115 * t97 - t149 * t186 - (qJD(3) * t130 * t137 + t135 * t154) * t116 * t136) * t110) * t105, 0, t148 * t168 * t181 + (t148 * t131 * t171 + (-t148 * t175 + ((-qJD(4) * t109 - 0.2e1 * t166) * t137 + (-t137 * t97 + (-qJD(4) * t116 + t98) * t135) * t110) * t138) * t130) * t105, t168 + 0.2e1 * (-t105 * t110 * t97 + (-t105 * t194 - t110 * t195) * t116) * t116, 0, 0;];
JaD_rot  = t1;
