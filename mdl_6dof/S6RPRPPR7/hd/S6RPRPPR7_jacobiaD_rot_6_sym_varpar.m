% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RPRPPR7
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR7_jacobiaD_rot_6_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_jacobiaD_rot_6_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_jacobiaD_rot_6_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:42:33
% EndTime: 2019-02-26 20:42:33
% DurationCPUTime: 0.69s
% Computational Cost: add. (1892->92), mult. (2519->206), div. (480->12), fcn. (2968->9), ass. (0->94)
t140 = cos(qJ(1));
t136 = t140 ^ 2;
t134 = qJ(3) + pkin(9);
t132 = sin(t134);
t128 = t132 ^ 2;
t133 = cos(t134);
t130 = 0.1e1 / t133 ^ 2;
t185 = t128 * t130;
t124 = t136 * t185 + 0.1e1;
t127 = t132 * t128;
t129 = 0.1e1 / t133;
t183 = t129 * t132;
t148 = qJD(3) * (t127 * t129 * t130 + t183);
t138 = sin(qJ(1));
t175 = qJD(1) * t140;
t152 = t128 * t138 * t175;
t192 = (-t130 * t152 + t136 * t148) / t124 ^ 2;
t203 = -0.2e1 * t192;
t162 = 0.1e1 + t185;
t202 = t140 * t162;
t155 = qJD(1) * t133 + qJD(6);
t172 = qJD(3) * t140;
t201 = t132 * t172 + t155 * t138;
t177 = t140 * t132;
t123 = atan2(t177, t133);
t119 = sin(t123);
t120 = cos(t123);
t105 = t119 * t177 + t120 * t133;
t102 = 0.1e1 / t105;
t139 = cos(qJ(6));
t178 = t139 * t140;
t137 = sin(qJ(6));
t180 = t138 * t137;
t118 = -t133 * t180 + t178;
t112 = 0.1e1 / t118;
t103 = 0.1e1 / t105 ^ 2;
t113 = 0.1e1 / t118 ^ 2;
t121 = 0.1e1 / t124;
t200 = t121 - 0.1e1;
t135 = t138 ^ 2;
t174 = qJD(3) * t133;
t184 = t128 * t135;
t164 = t130 * t172;
t176 = qJD(1) * t138;
t165 = t132 * t176;
t96 = ((t133 * t172 - t165) * t129 + t128 * t164) * t121;
t157 = -t96 + t172;
t158 = t140 * t96 - qJD(3);
t187 = t120 * t132;
t91 = t158 * t187 + (t157 * t133 - t165) * t119;
t197 = t102 * t103 * t91;
t99 = t103 * t184 + 0.1e1;
t199 = (-t184 * t197 + (t132 * t135 * t174 + t152) * t103) / t99 ^ 2;
t97 = 0.1e1 / t99;
t198 = t103 * t97;
t156 = qJD(6) * t133 + qJD(1);
t166 = t133 * t178;
t100 = -qJD(1) * t166 - qJD(6) * t178 + (qJD(3) * t132 * t139 + t156 * t137) * t138;
t179 = t138 * t139;
t181 = t137 * t140;
t117 = t133 * t179 + t181;
t111 = t117 ^ 2;
t110 = t111 * t113 + 0.1e1;
t190 = t113 * t117;
t173 = qJD(3) * t138;
t101 = -t156 * t179 + (t132 * t173 - t155 * t140) * t137;
t194 = t101 * t112 * t113;
t196 = 0.1e1 / t110 ^ 2 * (-t100 * t190 - t111 * t194);
t191 = t112 * t139;
t189 = t117 * t137;
t188 = t119 * t133;
t186 = t128 * t129;
t182 = t132 * t138;
t171 = 0.2e1 * t199;
t170 = 0.2e1 * t197;
t169 = 0.2e1 * t196;
t168 = t97 * t174;
t167 = t140 * t186;
t161 = -0.2e1 * t102 * t199;
t160 = 0.2e1 * t117 * t194;
t159 = 0.2e1 * t132 * t192;
t154 = t121 * t167;
t153 = t200 * t132 * t119;
t151 = t162 * t138;
t150 = t156 * t140;
t149 = t113 * t189 + t191;
t116 = -t133 * t181 - t179;
t115 = t166 - t180;
t108 = 0.1e1 / t110;
t107 = t121 * t202;
t95 = (-t120 * t154 + t153) * t138;
t93 = t140 * t188 - t187 + (t120 * t177 - t188) * t107;
t92 = t202 * t203 + (-qJD(1) * t151 + 0.2e1 * t140 * t148) * t121;
t1 = [t138 * t129 * t159 + (-qJD(3) * t151 - t175 * t183) * t121, 0, t92, 0, 0, 0; (t102 * t168 + (t161 + (-qJD(1) * t95 - t91) * t198) * t132) * t140 + ((-t95 * t168 + (t95 * t171 + ((-t96 * t154 - t200 * t174 + t159) * t119 + (t167 * t203 + t132 * t96 + (t127 * t164 - (t96 - 0.2e1 * t172) * t132) * t121) * t120) * t97 * t138) * t132) * t103 + (t95 * t170 + (-t102 + ((-t135 + t136) * t121 * t120 * t186 - t140 * t153) * t103) * qJD(1)) * t132 * t97) * t138, 0 (t102 * t97 * t175 + (t161 + (-qJD(3) * t93 - t91) * t198) * t138) * t133 + (((-qJD(3) * t102 + t93 * t170) * t138 + (-t93 * t175 + (-(-t107 * t176 + t140 * t92) * t120 - ((-t107 * t140 + 0.1e1) * t96 + (t107 - t140) * qJD(3)) * t119) * t182) * t103) * t97 + (t93 * t171 - ((-t92 - t176) * t119 + (t157 * t107 + t158) * t120) * t97 * t133) * t103 * t138) * t132, 0, 0, 0; (-t112 * t115 + t116 * t190) * t169 + (t116 * t160 - t112 * t137 * t150 - t201 * t191 + (t117 * t139 * t150 + t116 * t100 - t115 * t101 - t201 * t189) * t113) * t108, 0, t149 * t169 * t182 + (-t149 * t133 * t173 + (-t149 * t175 + ((qJD(6) * t112 + t160) * t137 + (t100 * t137 + (-qJD(6) * t117 + t101) * t139) * t113) * t138) * t132) * t108, 0, 0, -0.2e1 * t196 + 0.2e1 * (-t100 * t108 * t113 + (-t108 * t194 - t113 * t196) * t117) * t117;];
JaD_rot  = t1;
