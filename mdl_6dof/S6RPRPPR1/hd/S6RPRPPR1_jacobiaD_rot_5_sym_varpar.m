% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR1
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:39
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR1_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR1_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPPR1_jacobiaD_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:39:11
% EndTime: 2019-02-26 20:39:12
% DurationCPUTime: 0.69s
% Computational Cost: add. (2921->86), mult. (2191->195), div. (456->12), fcn. (2616->9), ass. (0->91)
t119 = qJ(1) + pkin(9);
t115 = sin(t119);
t109 = t115 ^ 2;
t118 = qJ(3) + pkin(10);
t114 = sin(t118);
t108 = t114 ^ 2;
t116 = cos(t118);
t111 = 0.1e1 / t116 ^ 2;
t161 = t108 * t111;
t106 = t109 * t161 + 0.1e1;
t107 = t114 * t108;
t110 = 0.1e1 / t116;
t160 = t110 * t114;
t129 = qJD(3) * (t107 * t110 * t111 + t160);
t117 = cos(t119);
t152 = qJD(1) * t117;
t139 = t115 * t152;
t166 = 0.1e1 / t106 ^ 2 * (t109 * t129 + t139 * t161);
t180 = -0.2e1 * t166;
t104 = 0.1e1 / t106;
t134 = 0.1e1 + t161;
t176 = t115 * t134;
t87 = t104 * t176;
t179 = t115 * t87 - 0.1e1;
t121 = cos(pkin(11));
t151 = qJD(3) * t114;
t136 = t121 * t151;
t120 = sin(pkin(11));
t155 = t117 * t120;
t156 = t115 * t121;
t98 = -t116 * t156 + t155;
t92 = qJD(1) * t98 - t117 * t136;
t154 = t117 * t121;
t157 = t115 * t120;
t100 = t116 * t154 + t157;
t95 = 0.1e1 / t100 ^ 2;
t178 = t92 * t95;
t99 = t116 * t155 - t156;
t168 = t95 * t99;
t93 = t99 ^ 2;
t90 = t93 * t95 + 0.1e1;
t88 = 0.1e1 / t90;
t94 = 0.1e1 / t100;
t177 = (-t120 * t94 + t121 * t168) * t88;
t158 = t115 * t114;
t103 = atan2(-t158, -t116);
t102 = cos(t103);
t101 = sin(t103);
t142 = t101 * t158;
t85 = -t102 * t116 - t142;
t82 = 0.1e1 / t85;
t83 = 0.1e1 / t85 ^ 2;
t175 = t104 - 0.1e1;
t113 = t117 ^ 2;
t149 = qJD(3) * t116;
t143 = t83 * t149;
t140 = t114 * t152;
t150 = qJD(3) * t115;
t162 = t102 * t114;
t138 = t111 * t150;
t78 = (-(-t115 * t149 - t140) * t110 + t108 * t138) * t104;
t73 = (-t115 * t78 + qJD(3)) * t162 + (-t140 + (t78 - t150) * t116) * t101;
t173 = t73 * t82 * t83;
t81 = t108 * t113 * t83 + 0.1e1;
t174 = (t113 * t114 * t143 + (-t113 * t173 - t139 * t83) * t108) / t81 ^ 2;
t169 = t94 * t178;
t137 = t120 * t151;
t97 = -t116 * t157 - t154;
t91 = qJD(1) * t97 - t117 * t137;
t172 = (t91 * t168 - t93 * t169) / t90 ^ 2;
t79 = 0.1e1 / t81;
t171 = t79 * t83;
t170 = t82 * t79;
t164 = t117 * t83;
t163 = qJD(3) * t87;
t159 = t114 * t117;
t153 = qJD(1) * t115;
t148 = qJD(3) * t117;
t147 = 0.2e1 * t173;
t146 = 0.2e1 * t172;
t145 = t82 * t174;
t144 = t99 * t169;
t141 = t104 * t108 * t110;
t135 = 0.2e1 * t83 * t174;
t133 = t110 * t180;
t132 = t115 * t141;
t131 = t134 * t117;
t77 = (t175 * t114 * t101 - t102 * t132) * t117;
t75 = -t179 * t162 + (-t115 + t87) * t116 * t101;
t74 = t176 * t180 + (qJD(1) * t131 + 0.2e1 * t115 * t129) * t104;
t1 = [t133 * t159 + (qJD(3) * t131 - t153 * t160) * t104, 0, t74, 0, 0, 0; (-t149 * t170 + (0.2e1 * t145 + (qJD(1) * t77 + t73) * t171) * t114) * t115 + (t77 * t135 * t114 + (-t77 * t143 + (t77 * t147 + ((0.2e1 * t114 * t166 - t78 * t132 - t175 * t149) * t101 + (t108 * t115 * t133 + t114 * t78 + (t107 * t138 - (t78 - 0.2e1 * t150) * t114) * t104) * t102) * t164) * t114 + (-t82 + (-(t109 - t113) * t102 * t141 + t175 * t142) * t83) * t114 * qJD(1)) * t79) * t117, 0 (-t153 * t170 + (-0.2e1 * t145 + (-qJD(3) * t75 - t73) * t171) * t117) * t116 + (t75 * t117 * t135 + (-t82 * t148 - ((-t115 * t74 - t152 * t87) * t102 + (t179 * t78 + t150 - t163) * t101) * t83 * t159 + (t117 * t147 + t83 * t153) * t75 - ((t74 - t152) * t101 + (t78 * t87 + qJD(3) + (-t78 - t163) * t115) * t102) * t116 * t164) * t79) * t114, 0, 0, 0; (t98 * t168 - t94 * t97) * t146 + ((-qJD(1) * t99 + t115 * t137) * t94 + 0.2e1 * t98 * t144 + (-t97 * t92 - (-qJD(1) * t100 + t115 * t136) * t99 - t98 * t91) * t95) * t88, 0, t116 * t148 * t177 + (-t153 * t177 + ((t94 * t146 + t88 * t178) * t120 + (-0.2e1 * t168 * t172 + (t91 * t95 - 0.2e1 * t144) * t88) * t121) * t117) * t114, 0, 0, 0;];
JaD_rot  = t1;
