% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPRP3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:44
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPRP3_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:44:42
% EndTime: 2019-02-26 20:44:42
% DurationCPUTime: 0.70s
% Computational Cost: add. (2324->95), mult. (2519->212), div. (480->12), fcn. (2968->9), ass. (0->95)
t132 = qJ(1) + pkin(9);
t128 = sin(t132);
t125 = t128 ^ 2;
t138 = sin(qJ(3));
t134 = t138 ^ 2;
t139 = cos(qJ(3));
t136 = 0.1e1 / t139 ^ 2;
t179 = t134 * t136;
t122 = t125 * t179 + 0.1e1;
t133 = t138 * t134;
t135 = 0.1e1 / t139;
t147 = qJD(3) * (t133 * t136 + t138) * t135;
t130 = cos(t132);
t177 = qJD(1) * t130;
t185 = t128 * t134;
t152 = t177 * t185;
t193 = (t125 * t147 + t136 * t152) / t122 ^ 2;
t202 = -0.2e1 * t193;
t131 = pkin(10) + qJ(5);
t127 = sin(t131);
t129 = cos(t131);
t180 = t130 * t139;
t115 = t128 * t127 + t129 * t180;
t110 = 0.1e1 / t115 ^ 2;
t186 = t128 * t129;
t114 = t127 * t180 - t186;
t190 = t114 * t129;
t109 = 0.1e1 / t115;
t192 = t109 * t127;
t149 = t110 * t190 - t192;
t108 = t114 ^ 2;
t101 = t108 * t110 + 0.1e1;
t99 = 0.1e1 / t101;
t201 = t149 * t99;
t159 = 0.1e1 + t179;
t200 = t128 * t159;
t184 = t128 * t138;
t119 = atan2(-t184, -t139);
t117 = cos(t119);
t116 = sin(t119);
t166 = t116 * t184;
t105 = -t117 * t139 - t166;
t102 = 0.1e1 / t105;
t103 = 0.1e1 / t105 ^ 2;
t120 = 0.1e1 / t122;
t199 = t120 - 0.1e1;
t126 = t130 ^ 2;
t173 = qJD(3) * t139;
t188 = t126 * t134;
t175 = qJD(3) * t128;
t161 = t136 * t175;
t176 = qJD(1) * t138;
t162 = t130 * t176;
t95 = (-(-t128 * t173 - t162) * t135 + t134 * t161) * t120;
t156 = t95 - t175;
t157 = -t128 * t95 + qJD(3);
t189 = t117 * t138;
t89 = t157 * t189 + (t139 * t156 - t162) * t116;
t195 = t102 * t103 * t89;
t98 = t103 * t188 + 0.1e1;
t198 = (-t188 * t195 + (t126 * t138 * t173 - t152) * t103) / t98 ^ 2;
t191 = t110 * t114;
t154 = -qJD(1) * t139 + qJD(5);
t155 = qJD(5) * t139 - qJD(1);
t174 = qJD(3) * t138;
t160 = t130 * t174;
t182 = t130 * t127;
t94 = -t155 * t182 + (t128 * t154 - t160) * t129;
t194 = t109 * t110 * t94;
t183 = t128 * t139;
t148 = t127 * t183 + t130 * t129;
t93 = t148 * qJD(1) - t115 * qJD(5) + t127 * t160;
t197 = 0.1e1 / t101 ^ 2 * (-t108 * t194 - t191 * t93);
t96 = 0.1e1 / t98;
t196 = t103 * t96;
t181 = t130 * t138;
t178 = qJD(1) * t128;
t172 = 0.2e1 * t198;
t171 = -0.2e1 * t197;
t170 = 0.2e1 * t195;
t169 = t102 * t198;
t168 = t114 * t194;
t167 = t96 * t173;
t165 = t120 * t134 * t135;
t163 = t128 * t176;
t158 = t135 * t202;
t153 = t128 * t165;
t151 = t159 * t130;
t150 = t154 * t130;
t113 = -t129 * t183 + t182;
t107 = t120 * t200;
t92 = (t116 * t138 * t199 - t117 * t153) * t130;
t91 = -t116 * t183 + t189 + (t116 * t139 - t117 * t184) * t107;
t90 = t200 * t202 + (qJD(1) * t151 + 0.2e1 * t128 * t147) * t120;
t1 = [t158 * t181 + (qJD(3) * t151 - t135 * t163) * t120, 0, t90, 0, 0, 0; (-t102 * t167 + (0.2e1 * t169 + (qJD(1) * t92 + t89) * t196) * t138) * t128 + (t92 * t138 * t96 * t170 + (-t92 * t167 + (t92 * t172 + ((0.2e1 * t138 * t193 - t153 * t95 - t173 * t199) * t116 + (t158 * t185 + t138 * t95 + (t133 * t161 - (t95 - 0.2e1 * t175) * t138) * t120) * t117) * t96 * t130) * t138) * t103 + (-t102 + ((-t125 + t126) * t117 * t165 + t199 * t166) * t103) * t96 * t176) * t130, 0 (-t102 * t96 * t178 + (-0.2e1 * t169 + (-qJD(3) * t91 - t89) * t196) * t130) * t139 + (t91 * t130 * t103 * t172 + ((-qJD(3) * t102 + t170 * t91) * t130 + (t91 * t178 + (-(-t107 * t177 - t128 * t90) * t117 - ((t107 * t128 - 0.1e1) * t95 + (-t107 + t128) * qJD(3)) * t116) * t181) * t103) * t96 - ((t90 - t177) * t116 + (t107 * t156 + t157) * t117) * t180 * t196) * t138, 0, 0, 0; 0.2e1 * (t109 * t148 + t113 * t191) * t197 + (0.2e1 * t113 * t168 - t155 * t109 * t186 + (t128 * t174 + t150) * t192 + (t148 * t94 + t113 * t93 - t150 * t190 - (t127 * t155 + t129 * t174) * t114 * t128) * t110) * t99, 0, -t163 * t201 + (t173 * t201 + (t149 * t171 + ((-qJD(5) * t109 - 0.2e1 * t168) * t129 + (-t129 * t93 + (-qJD(5) * t114 + t94) * t127) * t110) * t99) * t138) * t130, 0, t171 + 0.2e1 * (-t110 * t93 * t99 + (-t110 * t197 - t194 * t99) * t114) * t114, 0;];
JaD_rot  = t1;
