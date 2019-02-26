% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRRPP6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:59
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRRPP6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP6_jacobiaD_rot_5_sym_varpar: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:59:25
% EndTime: 2019-02-26 20:59:26
% DurationCPUTime: 0.68s
% Computational Cost: add. (1161->91), mult. (2519->207), div. (480->12), fcn. (2968->9), ass. (0->95)
t137 = cos(qJ(1));
t200 = 0.2e1 * t137;
t126 = qJ(4) + pkin(9);
t124 = sin(t126);
t134 = sin(qJ(3));
t125 = cos(t126);
t135 = sin(qJ(1));
t181 = t135 * t125;
t114 = t124 * t137 + t134 * t181;
t111 = 0.1e1 / t114 ^ 2;
t179 = t137 * t125;
t182 = t135 * t124;
t113 = t134 * t182 - t179;
t188 = t113 * t125;
t110 = 0.1e1 / t114;
t190 = t110 * t124;
t148 = t111 * t188 - t190;
t109 = t113 ^ 2;
t102 = t109 * t111 + 0.1e1;
t99 = 0.1e1 / t102;
t199 = t148 * t99;
t136 = cos(qJ(3));
t178 = t137 * t136;
t119 = atan2(-t178, t134);
t117 = sin(t119);
t118 = cos(t119);
t106 = -t117 * t178 + t118 * t134;
t103 = 0.1e1 / t106;
t127 = 0.1e1 / t134;
t104 = 0.1e1 / t106 ^ 2;
t128 = 0.1e1 / t134 ^ 2;
t133 = t137 ^ 2;
t132 = t136 ^ 2;
t185 = t128 * t132;
t122 = t133 * t185 + 0.1e1;
t120 = 0.1e1 / t122;
t198 = t120 - 0.1e1;
t130 = t135 ^ 2;
t184 = t130 * t132;
t101 = t104 * t184 + 0.1e1;
t175 = qJD(1) * t137;
t152 = t132 * t135 * t175;
t173 = qJD(3) * t136;
t172 = qJD(3) * t137;
t162 = t128 * t172;
t176 = qJD(1) * t136;
t164 = t135 * t176;
t96 = ((t134 * t172 + t164) * t127 + t132 * t162) * t120;
t157 = -t96 + t172;
t158 = -t137 * t96 + qJD(3);
t187 = t118 * t136;
t90 = t158 * t187 + (t157 * t134 + t164) * t117;
t194 = t103 * t104 * t90;
t197 = (-t184 * t194 + (-t130 * t134 * t173 + t152) * t104) / t101 ^ 2;
t189 = t111 * t113;
t155 = qJD(1) * t134 + qJD(4);
t146 = t135 * t173 + t155 * t137;
t156 = qJD(4) * t134 + qJD(1);
t150 = t124 * t156;
t95 = t146 * t125 - t135 * t150;
t193 = t110 * t111 * t95;
t149 = t125 * t156;
t94 = t146 * t124 + t135 * t149;
t196 = 0.1e1 / t102 ^ 2 * (-t109 * t193 + t94 * t189);
t97 = 0.1e1 / t101;
t195 = t104 * t97;
t131 = t136 * t132;
t147 = qJD(3) * (-t128 * t131 - t136) * t127;
t191 = (-t128 * t152 + t133 * t147) / t122 ^ 2;
t186 = t127 * t132;
t183 = t134 * t137;
t177 = qJD(1) * t135;
t174 = qJD(3) * t134;
t171 = -0.2e1 * t197;
t170 = 0.2e1 * t196;
t169 = -0.2e1 * t194;
t168 = t103 * t197;
t167 = t97 * t174;
t166 = t136 * t191;
t165 = t120 * t186;
t163 = t136 * t175;
t161 = 0.1e1 + t185;
t160 = 0.2e1 * t113 * t193;
t159 = t191 * t200;
t154 = t137 * t165;
t153 = t198 * t136 * t117;
t151 = t161 * t135;
t145 = -t155 * t135 + t136 * t172;
t116 = t134 * t179 - t182;
t115 = t124 * t183 + t181;
t108 = t161 * t137 * t120;
t93 = (-t118 * t154 - t153) * t135;
t92 = t117 * t183 + t187 + (-t117 * t134 - t118 * t178) * t108;
t91 = -t161 * t159 + (-qJD(1) * t151 + t147 * t200) * t120;
t1 = [-0.2e1 * t135 * t127 * t166 + (-qJD(3) * t151 + t127 * t163) * t120, 0, t91, 0, 0, 0; (t103 * t167 + (0.2e1 * t168 + (qJD(1) * t93 + t90) * t195) * t136) * t137 + (t93 * t136 * t97 * t169 + (-t93 * t167 + (t93 * t171 + ((t96 * t154 + t198 * t174 + 0.2e1 * t166) * t117 + (t159 * t186 + t96 * t136 + (t131 * t162 + (-t96 + 0.2e1 * t172) * t136) * t120) * t118) * t97 * t135) * t136) * t104 + (t103 + ((t130 - t133) * t118 * t165 - t137 * t153) * t104) * t97 * t176) * t135, 0 (t103 * t97 * t175 + (-0.2e1 * t168 + (-qJD(3) * t92 - t90) * t195) * t135) * t134 + (((qJD(3) * t103 + t92 * t169) * t135 + (t92 * t175 + ((t108 * t177 - t137 * t91) * t118 + ((t108 * t137 - 0.1e1) * t96 + (-t108 + t137) * qJD(3)) * t117) * t135 * t136) * t104) * t97 + (t92 * t171 + ((-t91 - t177) * t117 + (t157 * t108 - t158) * t118) * t97 * t134) * t104 * t135) * t136, 0, 0, 0; (-t110 * t115 + t116 * t189) * t170 + (t116 * t160 + t137 * t110 * t149 + t145 * t190 + (t137 * t113 * t150 - t115 * t95 - t116 * t94 - t145 * t188) * t111) * t99, 0, -t163 * t199 + (t174 * t199 + (t148 * t170 + ((qJD(4) * t110 + t160) * t125 + (-t125 * t94 + (qJD(4) * t113 - t95) * t124) * t111) * t99) * t136) * t135, -0.2e1 * t196 + 0.2e1 * (t111 * t94 * t99 + (-t111 * t196 - t99 * t193) * t113) * t113, 0, 0;];
JaD_rot  = t1;
