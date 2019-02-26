% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RPRPPR6
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4,theta5]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:42
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RPRPPR6_jacobiaD_rot_5_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR6_jacobiaD_rot_5_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR6_jacobiaD_rot_5_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:41:58
% EndTime: 2019-02-26 20:41:59
% DurationCPUTime: 0.67s
% Computational Cost: add. (1897->83), mult. (2191->188), div. (456->12), fcn. (2616->9), ass. (0->86)
t107 = qJ(3) + pkin(9);
t106 = cos(t107);
t104 = t106 ^ 2;
t105 = sin(t107);
t150 = 0.1e1 / t105 ^ 2 * t104;
t113 = cos(qJ(1));
t126 = 0.1e1 + t150;
t109 = t113 ^ 2;
t99 = t109 * t150 + 0.1e1;
t97 = 0.1e1 / t99;
t153 = t113 * t97;
t80 = t126 * t153;
t170 = t113 * t80 - 0.1e1;
t111 = cos(pkin(10));
t112 = sin(qJ(1));
t140 = qJD(3) * t112;
t130 = t106 * t140;
t110 = sin(pkin(10));
t147 = t112 * t110;
t148 = t111 * t113;
t93 = t105 * t148 - t147;
t85 = t93 * qJD(1) + t111 * t130;
t146 = t112 * t111;
t149 = t110 * t113;
t91 = t105 * t146 + t149;
t88 = 0.1e1 / t91 ^ 2;
t169 = t85 * t88;
t90 = t105 * t147 - t148;
t157 = t88 * t90;
t86 = t90 ^ 2;
t83 = t86 * t88 + 0.1e1;
t81 = 0.1e1 / t83;
t87 = 0.1e1 / t91;
t168 = (-t110 * t87 + t111 * t157) * t81;
t167 = t106 * t150;
t145 = t113 * t106;
t96 = atan2(-t145, t105);
t94 = sin(t96);
t95 = cos(t96);
t78 = t95 * t105 - t94 * t145;
t75 = 0.1e1 / t78;
t100 = 0.1e1 / t105;
t76 = 0.1e1 / t78 ^ 2;
t166 = t97 - 0.1e1;
t108 = t112 ^ 2;
t142 = qJD(1) * t113;
t131 = t112 * t142;
t141 = qJD(3) * t105;
t132 = t76 * t141;
t139 = qJD(3) * t113;
t143 = qJD(1) * t112;
t156 = t105 * t94;
t71 = ((t105 * t139 + t106 * t143) * t100 + t139 * t150) * t97;
t66 = (-t71 + t139) * t156 + (t94 * t143 + (-t113 * t71 + qJD(3)) * t95) * t106;
t164 = t66 * t75 * t76;
t74 = t104 * t108 * t76 + 0.1e1;
t165 = (-t108 * t106 * t132 + (-t108 * t164 + t76 * t131) * t104) / t74 ^ 2;
t158 = t87 * t169;
t92 = t105 * t149 + t146;
t84 = t92 * qJD(1) + t110 * t130;
t163 = (t84 * t157 - t86 * t158) / t83 ^ 2;
t72 = 0.1e1 / t74;
t161 = t72 * t76;
t160 = t75 * t72;
t121 = qJD(3) * (-t106 - t167) * t100;
t159 = (t109 * t121 - t131 * t150) / t99 ^ 2;
t155 = t112 * t76;
t152 = qJD(3) * t80;
t151 = t100 * t104;
t144 = qJD(1) * t106;
t138 = -0.2e1 * t164;
t137 = 0.2e1 * t163;
t136 = t75 * t165;
t135 = t106 * t159;
t134 = t100 * t153;
t133 = t106 * t166;
t129 = t106 * t139;
t128 = -0.2e1 * t76 * t165;
t127 = 0.2e1 * t90 * t158;
t125 = t94 * t133;
t124 = t104 * t134;
t123 = t126 * t97;
t70 = (-t124 * t95 - t125) * t112;
t68 = -t170 * t95 * t106 + (t113 - t80) * t156;
t67 = -t123 * t143 + 0.2e1 * (t121 * t97 - t126 * t159) * t113;
t1 = [t134 * t144 + (-qJD(3) * t123 - 0.2e1 * t100 * t135) * t112, 0, t67, 0, 0, 0; (t141 * t160 + (0.2e1 * t136 + (qJD(1) * t70 + t66) * t161) * t106) * t113 + (t70 * t128 * t106 + (-t70 * t132 + (t70 * t138 + ((t71 * t124 + t166 * t141 + 0.2e1 * t135) * t94 + (-t71 * t133 + (0.2e1 * t151 * t159 + (0.2e1 * t106 + t167) * t97 * qJD(3)) * t113) * t95) * t155) * t106 + (t75 + (-t113 * t125 + (t108 - t109) * t97 * t95 * t151) * t76) * t144) * t72) * t112, 0 (t142 * t160 + (-0.2e1 * t136 + (-qJD(3) * t68 - t66) * t161) * t112) * t105 + (t68 * t112 * t128 + (t75 * t140 + (t112 * t138 + t76 * t142) * t68 + (((-t113 * t67 + t143 * t80) * t95 + (t170 * t71 + t139 - t152) * t94) * t106 + ((-t67 - t143) * t94 + (-t71 * t80 - qJD(3) + (t71 + t152) * t113) * t95) * t105) * t155) * t72) * t106, 0, 0, 0; (t93 * t157 - t87 * t92) * t137 + ((-t90 * qJD(1) + t110 * t129) * t87 + t93 * t127 + (-t92 * t85 - (-t91 * qJD(1) + t111 * t129) * t90 - t93 * t84) * t88) * t81, 0, t105 * t140 * t168 + (-t142 * t168 + ((-0.2e1 * t87 * t163 - t81 * t169) * t110 + (t137 * t157 + (-t84 * t88 + t127) * t81) * t111) * t112) * t106, 0, 0, 0;];
JaD_rot  = t1;
