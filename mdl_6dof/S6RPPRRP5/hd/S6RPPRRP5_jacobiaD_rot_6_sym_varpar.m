% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:32
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPPRRP5_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_jacobiaD_rot_6_sym_varpar: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:32:33
% EndTime: 2019-02-26 20:32:33
% DurationCPUTime: 0.66s
% Computational Cost: add. (624->92), mult. (2519->208), div. (480->12), fcn. (2968->9), ass. (0->94)
t134 = sin(qJ(1));
t133 = sin(qJ(4));
t126 = 0.1e1 / t133 ^ 2;
t136 = cos(qJ(4));
t130 = t136 ^ 2;
t186 = t126 * t130;
t158 = 0.1e1 + t186;
t198 = t134 * t158;
t179 = t134 * t136;
t120 = atan2(t179, t133);
t117 = cos(t120);
t116 = sin(t120);
t165 = t116 * t179;
t102 = t117 * t133 + t165;
t99 = 0.1e1 / t102;
t135 = cos(qJ(5));
t137 = cos(qJ(1));
t178 = t135 * t137;
t163 = t133 * t178;
t132 = sin(qJ(5));
t181 = t134 * t132;
t115 = t163 - t181;
t109 = 0.1e1 / t115;
t125 = 0.1e1 / t133;
t100 = 0.1e1 / t102 ^ 2;
t110 = 0.1e1 / t115 ^ 2;
t128 = t134 ^ 2;
t121 = t128 * t186 + 0.1e1;
t118 = 0.1e1 / t121;
t197 = t118 - 0.1e1;
t131 = t137 ^ 2;
t185 = t130 * t131;
t98 = t100 * t185 + 0.1e1;
t96 = 0.1e1 / t98;
t196 = t100 * t96;
t172 = qJD(4) * t134;
t161 = t126 * t172;
t174 = qJD(1) * t137;
t162 = t136 * t174;
t93 = ((-t133 * t172 + t162) * t125 - t130 * t161) * t118;
t154 = -t93 - t172;
t155 = t134 * t93 + qJD(4);
t187 = t117 * t136;
t88 = t155 * t187 + (t154 * t133 + t162) * t116;
t195 = t99 * t100 * t88;
t180 = t134 * t135;
t182 = t133 * t137;
t114 = t132 * t182 + t180;
t108 = t114 ^ 2;
t107 = t108 * t110 + 0.1e1;
t190 = t110 * t114;
t152 = qJD(1) * t133 + qJD(5);
t153 = qJD(5) * t133 + qJD(1);
t171 = qJD(4) * t136;
t160 = t137 * t171;
t183 = t132 * t137;
t95 = -t153 * t183 + (-t152 * t134 + t160) * t135;
t193 = t109 * t110 * t95;
t94 = -qJD(5) * t163 - t132 * t160 - t135 * t174 + t152 * t181;
t194 = 0.1e1 / t107 ^ 2 * (-t108 * t193 - t94 * t190);
t129 = t136 * t130;
t148 = (t126 * t129 + t136) * t125;
t184 = t130 * t134;
t150 = t174 * t184;
t192 = (-t148 * t128 * qJD(4) + t126 * t150) / t121 ^ 2;
t191 = t109 * t132;
t189 = t114 * t135;
t188 = t116 * t133;
t177 = t136 * t137;
t176 = qJD(1) * t134;
t175 = qJD(1) * t136;
t173 = qJD(4) * t133;
t170 = -0.2e1 * (-t185 * t195 + (-t131 * t133 * t171 - t150) * t100) / t98 ^ 2;
t169 = -0.2e1 * t195;
t168 = 0.2e1 * t194;
t167 = 0.2e1 * t192;
t166 = t96 * t173;
t164 = t118 * t125 * t130;
t159 = t99 * t170;
t157 = 0.2e1 * t114 * t193;
t156 = -0.2e1 * t125 * t192;
t151 = t134 * t164;
t149 = t158 * t137;
t147 = t152 * t137;
t146 = t110 * t189 - t191;
t145 = t146 * t137;
t113 = -t133 * t180 - t183;
t112 = -t133 * t181 + t178;
t106 = t118 * t198;
t104 = 0.1e1 / t107;
t92 = (-t197 * t136 * t116 + t117 * t151) * t137;
t91 = -t134 * t188 + t187 - (t117 * t179 - t188) * t106;
t89 = t167 * t198 + (-qJD(1) * t149 + 0.2e1 * t148 * t172) * t118;
t1 = [t156 * t177 + (-t125 * t134 * t175 - qJD(4) * t149) * t118, 0, 0, t89, 0, 0; (-t99 * t166 + (t159 + (-qJD(1) * t92 - t88) * t196) * t136) * t134 + (t92 * t136 * t96 * t169 + (-t92 * t166 + (t92 * t170 + ((t136 * t167 - t93 * t151 + t197 * t173) * t116 + (t156 * t184 + t136 * t93 + (-t129 * t161 + (-t93 - 0.2e1 * t172) * t136) * t118) * t117) * t96 * t137) * t136) * t100 + (t99 + ((-t128 + t131) * t117 * t164 + t197 * t165) * t100) * t96 * t175) * t137, 0, 0 (-t99 * t96 * t176 + (t159 + (-qJD(4) * t91 - t88) * t196) * t137) * t133 + (t91 * t137 * t100 * t170 + ((qJD(4) * t99 + t91 * t169) * t137 + (-t91 * t176 + ((-t106 * t174 + t134 * t89) * t117 + ((t106 * t134 - 0.1e1) * t93 + (t106 - t134) * qJD(4)) * t116) * t177) * t100) * t96 + ((-t89 - t174) * t116 + (-t154 * t106 - t155) * t117) * t182 * t196) * t136, 0, 0; (-t109 * t112 + t113 * t190) * t168 + (t113 * t157 - t153 * t109 * t180 + (-t134 * t171 - t147) * t191 + (-t112 * t95 + t113 * t94 + t147 * t189 - (t153 * t132 - t135 * t171) * t114 * t134) * t110) * t104, 0, 0, t136 * t145 * t168 + (t145 * t173 + (t146 * t176 + ((qJD(5) * t109 + t157) * t135 + (t135 * t94 + (qJD(5) * t114 - t95) * t132) * t110) * t137) * t136) * t104, -0.2e1 * t194 + 0.2e1 * (-t104 * t110 * t94 + (-t104 * t193 - t110 * t194) * t114) * t114, 0;];
JaD_rot  = t1;
