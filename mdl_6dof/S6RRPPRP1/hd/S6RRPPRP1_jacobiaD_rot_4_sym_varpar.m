% Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6RRPPRP1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d5,theta3,theta4]';
%
% Output:
% JaD_rot [3x6]
%   Zeitableitung der rotatorischen Teilmatrix der analytischen Jacobi-Matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 21:24
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JaD_rot = S6RRPPRP1_jacobiaD_rot_4_sym_varpar(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPRP1_jacobiaD_rot_4_sym_varpar: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPRP1_jacobiaD_rot_4_sym_varpar: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From jacobiaD_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 21:24:45
% EndTime: 2019-02-26 21:24:46
% DurationCPUTime: 0.67s
% Computational Cost: add. (2086->84), mult. (2191->187), div. (456->12), fcn. (2616->9), ass. (0->87)
t107 = qJ(2) + pkin(9);
t105 = sin(t107);
t101 = t105 ^ 2;
t106 = cos(t107);
t150 = t101 / t106 ^ 2;
t112 = sin(qJ(1));
t126 = 0.1e1 + t150;
t108 = t112 ^ 2;
t99 = t108 * t150 + 0.1e1;
t97 = 0.1e1 / t99;
t123 = t126 * t97;
t80 = t112 * t123;
t171 = t112 * t80 - 0.1e1;
t111 = cos(pkin(10));
t113 = cos(qJ(1));
t139 = qJD(2) * t113;
t128 = t105 * t139;
t145 = t112 * t111;
t110 = sin(pkin(10));
t149 = t110 * t113;
t91 = -t106 * t145 + t149;
t85 = t91 * qJD(1) - t111 * t128;
t146 = t112 * t110;
t148 = t111 * t113;
t93 = t106 * t148 + t146;
t88 = 0.1e1 / t93 ^ 2;
t170 = t85 * t88;
t92 = t106 * t149 - t145;
t157 = t88 * t92;
t86 = t92 ^ 2;
t83 = t86 * t88 + 0.1e1;
t81 = 0.1e1 / t83;
t87 = 0.1e1 / t93;
t169 = (-t110 * t87 + t111 * t157) * t81;
t168 = t105 * t150;
t147 = t112 * t105;
t96 = atan2(-t147, -t106);
t94 = sin(t96);
t133 = t94 * t147;
t95 = cos(t96);
t78 = -t106 * t95 - t133;
t75 = 0.1e1 / t78;
t102 = 0.1e1 / t106;
t76 = 0.1e1 / t78 ^ 2;
t167 = 0.2e1 * t105;
t166 = t97 - 0.1e1;
t109 = t113 ^ 2;
t142 = qJD(1) * t113;
t130 = t112 * t142;
t141 = qJD(2) * t106;
t131 = t76 * t141;
t140 = qJD(2) * t112;
t154 = t106 * t94;
t71 = (-(-t105 * t142 - t106 * t140) * t102 + t140 * t150) * t97;
t66 = (t71 - t140) * t154 + (-t94 * t142 + (-t112 * t71 + qJD(2)) * t95) * t105;
t164 = t66 * t75 * t76;
t156 = t101 * t76;
t74 = t109 * t156 + 0.1e1;
t165 = (t109 * t105 * t131 + (-t109 * t164 - t76 * t130) * t101) / t74 ^ 2;
t158 = t87 * t170;
t90 = -t106 * t146 - t148;
t84 = t90 * qJD(1) - t110 * t128;
t163 = (t84 * t157 - t86 * t158) / t83 ^ 2;
t72 = 0.1e1 / t74;
t161 = t72 * t76;
t160 = t75 * t72;
t121 = qJD(2) * (t105 + t168) * t102;
t159 = (t108 * t121 + t130 * t150) / t99 ^ 2;
t155 = t102 * t97;
t152 = t113 * t76;
t151 = qJD(2) * t80;
t144 = qJD(1) * t105;
t143 = qJD(1) * t112;
t138 = 0.2e1 * t164;
t137 = 0.2e1 * t163;
t136 = t75 * t165;
t135 = t92 * t158;
t134 = t112 * t155;
t132 = t166 * t105;
t129 = t105 * t140;
t127 = 0.2e1 * t76 * t165;
t125 = -0.2e1 * t102 * t159;
t124 = t101 * t134;
t70 = (-t95 * t124 + t94 * t132) * t113;
t68 = (-t112 + t80) * t154 - t171 * t95 * t105;
t67 = t123 * t142 + 0.2e1 * (t121 * t97 - t126 * t159) * t112;
t1 = [-t134 * t144 + (qJD(2) * t123 + t105 * t125) * t113, t67, 0, 0, 0, 0; (-t141 * t160 + (0.2e1 * t136 + (qJD(1) * t70 + t66) * t161) * t105) * t112 + (t70 * t127 * t105 + (-t70 * t131 + (t70 * t138 + ((-t71 * t124 - t166 * t141 + t159 * t167) * t94 + (-t71 * t132 + (t101 * t125 + (t167 + t168) * t97 * qJD(2)) * t112) * t95) * t152) * t105 + (-t75 + t166 * t76 * t133 - (t108 - t109) * t95 * t155 * t156) * t144) * t72) * t113 (-t143 * t160 + (-0.2e1 * t136 + (-qJD(2) * t68 - t66) * t161) * t113) * t106 + (t68 * t113 * t127 + (-t75 * t139 + (t113 * t138 + t76 * t143) * t68 + (-((-t112 * t67 - t142 * t80) * t95 + (t171 * t71 + t140 - t151) * t94) * t105 - ((t67 - t142) * t94 + (t71 * t80 + qJD(2) + (-t71 - t151) * t112) * t95) * t106) * t152) * t72) * t105, 0, 0, 0, 0; (t91 * t157 - t87 * t90) * t137 + ((-t92 * qJD(1) + t110 * t129) * t87 + 0.2e1 * t91 * t135 + (-t90 * t85 - (-t93 * qJD(1) + t111 * t129) * t92 - t91 * t84) * t88) * t81, t106 * t139 * t169 + (-t143 * t169 + ((t87 * t137 + t81 * t170) * t110 + (-0.2e1 * t157 * t163 + (t84 * t88 - 0.2e1 * t135) * t81) * t111) * t113) * t105, 0, 0, 0, 0;];
JaD_rot  = t1;
