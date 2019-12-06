% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:10:13
% EndTime: 2019-12-05 17:10:18
% DurationCPUTime: 1.17s
% Computational Cost: add. (5007->188), mult. (6526->267), div. (0->0), fcn. (4395->10), ass. (0->135)
t135 = sin(qJ(5));
t129 = qJDD(4) + qJDD(5);
t132 = qJD(2) + qJD(3);
t139 = cos(qJ(5));
t140 = cos(qJ(4));
t136 = sin(qJ(4));
t169 = t132 * t136;
t97 = -t139 * t140 * t132 + t135 * t169;
t99 = (t140 * t135 + t136 * t139) * t132;
t75 = t99 * t97;
t182 = -t75 + t129;
t185 = t135 * t182;
t184 = t139 * t182;
t172 = sin(pkin(9));
t173 = cos(pkin(9));
t117 = -t172 * g(1) + t173 * g(2);
t128 = t132 ^ 2;
t130 = qJDD(2) + qJDD(3);
t137 = sin(qJ(3));
t141 = cos(qJ(3));
t138 = sin(qJ(2));
t142 = cos(qJ(2));
t148 = -t173 * g(1) - t172 * g(2);
t165 = -g(3) + qJDD(1);
t93 = -t138 * t148 + t142 * t165;
t146 = qJDD(2) * pkin(2) + t93;
t144 = qJD(2) ^ 2;
t94 = t138 * t165 + t142 * t148;
t92 = -t144 * pkin(2) + t94;
t67 = t137 * t146 + t141 * t92;
t61 = -t128 * pkin(3) + t130 * pkin(7) + t67;
t53 = -t140 * t117 + t136 * t61;
t54 = t136 * t117 + t140 * t61;
t30 = t136 * t53 + t140 * t54;
t164 = qJD(4) * t132;
t159 = t140 * t164;
t167 = t136 * t130;
t106 = t159 + t167;
t122 = t140 * t130;
t160 = t136 * t164;
t152 = t122 - t160;
t64 = -t97 * qJD(5) + t139 * t106 + t135 * t152;
t131 = qJD(4) + qJD(5);
t91 = t131 * t97;
t183 = -t91 + t64;
t181 = -t53 + (-t106 + t159) * pkin(8);
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t127 = t131 ^ 2;
t118 = t136 * t128 * t140;
t163 = qJDD(4) + t118;
t145 = t163 * pkin(4) + t181;
t116 = qJD(4) * pkin(4) - pkin(8) * t169;
t134 = t140 ^ 2;
t124 = t134 * t128;
t36 = -pkin(4) * t124 + t152 * pkin(8) - qJD(4) * t116 + t54;
t18 = t135 * t36 - t139 * t145;
t176 = t139 * t36;
t19 = t135 * t145 + t176;
t7 = t135 * t19 - t139 * t18;
t180 = t136 * t7;
t157 = t137 * t92 - t141 * t146;
t60 = -t130 * pkin(3) - t128 * pkin(7) + t157;
t179 = -pkin(3) * t60 + pkin(7) * t30;
t40 = -t152 * pkin(4) - pkin(8) * t124 + t116 * t169 + t60;
t178 = t135 * t40;
t72 = t75 + t129;
t177 = t135 * t72;
t175 = t139 * t40;
t174 = t139 * t72;
t171 = t131 * t135;
t170 = t131 * t139;
t168 = t136 * t163;
t166 = t140 * (qJDD(4) - t118);
t105 = 0.2e1 * t159 + t167;
t133 = t136 ^ 2;
t123 = t133 * t128;
t143 = qJD(4) ^ 2;
t85 = -t166 - t136 * (-t123 - t143);
t162 = -pkin(3) * t105 + pkin(7) * t85 + t136 * t60;
t107 = t122 - 0.2e1 * t160;
t84 = t140 * (-t124 - t143) - t168;
t161 = pkin(3) * t107 + pkin(7) * t84 - t140 * t60;
t154 = t135 * t106 - t139 * t152;
t147 = (-qJD(5) + t131) * t99 - t154;
t50 = t91 + t64;
t27 = t135 * t147 - t139 * t50;
t28 = t135 * t50 + t139 * t147;
t12 = -t136 * t27 + t140 * t28;
t65 = -t95 - t96;
t8 = t135 * t18 + t139 * t19;
t158 = t136 * (-pkin(8) * t27 - t7) + t140 * (-pkin(4) * t65 + pkin(8) * t28 + t8) - pkin(3) * t65 + pkin(7) * t12;
t70 = -t127 - t95;
t41 = t135 * t70 + t184;
t42 = t139 * t70 - t185;
t23 = -t136 * t41 + t140 * t42;
t45 = (qJD(5) + t131) * t99 + t154;
t156 = t136 * (-pkin(8) * t41 + t178) + t140 * (-pkin(4) * t45 + pkin(8) * t42 - t175) - pkin(3) * t45 + pkin(7) * t23;
t86 = -t96 - t127;
t55 = t139 * t86 - t177;
t56 = -t135 * t86 - t174;
t32 = -t136 * t55 + t140 * t56;
t155 = t136 * (-pkin(8) * t55 + t175) + t140 * (-pkin(4) * t183 + pkin(8) * t56 + t178) - pkin(3) * t183 + pkin(7) * t32;
t109 = (t133 + t134) * t130;
t112 = t123 + t124;
t153 = pkin(3) * t112 + pkin(7) * t109 + t30;
t110 = -t141 * t128 - t137 * t130;
t150 = t137 * t128 - t141 * t130;
t3 = t140 * t8 - t180;
t149 = pkin(7) * t3 - pkin(8) * t180 - pkin(3) * t40 + t140 * (-pkin(4) * t40 + pkin(8) * t8);
t89 = -t96 + t127;
t88 = t95 - t127;
t83 = t168 + t140 * (-t123 + t143);
t82 = t136 * (t124 - t143) + t166;
t79 = (t106 + t159) * t136;
t78 = t107 * t140;
t77 = t137 * t109 + t141 * t112;
t76 = t140 * t105 + t136 * t107;
t74 = t96 - t95;
t69 = -t141 * t105 + t137 * t85;
t68 = t141 * t107 + t137 * t84;
t63 = -t99 * qJD(5) - t154;
t38 = (t136 * (t135 * t99 - t139 * t97) + t140 * (-t135 * t97 - t139 * t99)) * t131;
t37 = t137 * t67 - t141 * t157;
t34 = t136 * (t139 * t88 - t177) + t140 * (t135 * t88 + t174);
t33 = t136 * (-t135 * t89 + t184) + t140 * (t139 * t89 + t185);
t26 = t136 * (t139 * t64 - t99 * t171) + t140 * (t135 * t64 + t99 * t170);
t25 = t136 * (-t135 * t63 + t97 * t170) + t140 * (t139 * t63 + t97 * t171);
t20 = t137 * t30 - t141 * t60;
t16 = t137 * t32 - t141 * t183;
t14 = t137 * t23 - t141 * t45;
t11 = t136 * (-t135 * t183 - t139 * t45) + t140 * (-t135 * t45 + t139 * t183);
t9 = t137 * t12 - t141 * t65;
t1 = t137 * t3 - t141 * t40;
t2 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t165, 0, 0, 0, 0, 0, 0, t142 * qJDD(2) - t138 * t144, -t138 * qJDD(2) - t142 * t144, 0, t138 * t94 + t142 * t93, 0, 0, 0, 0, 0, 0, t138 * t110 - t142 * t150, t142 * t110 + t138 * t150, 0, t138 * (t137 * t157 + t141 * t67) + t142 * t37, 0, 0, 0, 0, 0, 0, t138 * (-t137 * t107 + t141 * t84) + t142 * t68, t138 * (t137 * t105 + t141 * t85) + t142 * t69, t138 * (t141 * t109 - t137 * t112) + t142 * t77, t138 * (t137 * t60 + t141 * t30) + t142 * t20, 0, 0, 0, 0, 0, 0, t138 * (t137 * t45 + t141 * t23) + t142 * t14, t138 * (t137 * t183 + t141 * t32) + t142 * t16, t138 * (t141 * t12 + t137 * t65) + t142 * t9, t138 * (t137 * t40 + t141 * t3) + t142 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t93, -t94, 0, 0, 0, 0, 0, 0, 0, t130, -pkin(2) * t150 - t157, pkin(2) * t110 - t67, 0, pkin(2) * t37, t79, t76, t83, t78, t82, 0, pkin(2) * t68 + t161, pkin(2) * t69 + t162, pkin(2) * t77 + t153, pkin(2) * t20 + t179, t26, t11, t33, t25, t34, t38, pkin(2) * t14 + t156, pkin(2) * t16 + t155, pkin(2) * t9 + t158, pkin(2) * t1 + t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, -t157, -t67, 0, 0, t79, t76, t83, t78, t82, 0, t161, t162, t153, t179, t26, t11, t33, t25, t34, t38, t156, t155, t158, t149; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t118, t123 - t124, t167, t118, t122, qJDD(4), -t53, -t54, 0, 0, t75, t74, t50, -t75, t147, t129, pkin(4) * t41 - t18, -t176 - t135 * t181 + (-t135 * t163 + t55) * pkin(4), pkin(4) * t27, pkin(4) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t50, -t75, t147, t129, -t18, -t19, 0, 0;];
tauJ_reg = t2;
