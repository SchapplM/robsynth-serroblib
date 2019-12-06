% Calculate inertial parameters regressor of inverse dynamics joint torque vector for
% S5PRPRP3
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% tau_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:34
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5PRPRP3_invdynJ_fixb_reg2_slag_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP3_invdynJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:33:47
% EndTime: 2019-12-05 15:33:50
% DurationCPUTime: 1.10s
% Computational Cost: add. (1085->220), mult. (2244->278), div. (0->0), fcn. (1587->10), ass. (0->133)
t130 = qJD(1) * qJD(2);
t90 = sin(qJ(2));
t92 = cos(qJ(2));
t168 = t90 * qJDD(1) + t92 * t130;
t140 = qJDD(2) * pkin(2);
t77 = t92 * qJDD(1);
t41 = -t90 * t130 + t140 + t77;
t84 = sin(pkin(8));
t86 = cos(pkin(8));
t16 = t168 * t86 + t84 * t41;
t12 = qJDD(2) * pkin(6) + t16;
t169 = qJD(3) * qJD(4) + t12;
t85 = sin(pkin(7));
t87 = cos(pkin(7));
t114 = g(1) * t87 + g(2) * t85;
t81 = qJ(2) + pkin(8);
t71 = sin(t81);
t109 = t114 * t71;
t72 = cos(t81);
t160 = g(3) * t72;
t167 = t160 - t109;
t144 = qJD(1) * t90;
t135 = t92 * qJD(1);
t60 = qJD(2) * pkin(2) + t135;
t30 = t86 * t144 + t84 * t60;
t22 = qJD(2) * pkin(6) + t30;
t116 = qJ(5) * qJD(2) + t22;
t91 = cos(qJ(4));
t111 = t116 * t91;
t125 = -g(1) * t85 + g(2) * t87;
t161 = g(3) * t71;
t101 = t114 * t72 + t161;
t166 = pkin(2) * t90;
t89 = sin(qJ(4));
t82 = t89 ^ 2;
t165 = pkin(4) * t82;
t159 = g(3) * t89;
t158 = t86 * pkin(2);
t157 = t91 * pkin(4);
t156 = t85 * t89;
t155 = t85 * t91;
t154 = t87 * t89;
t153 = t87 * t91;
t63 = t84 * t144;
t40 = t86 * t135 - t63;
t152 = t91 * t40;
t94 = qJD(2) ^ 2;
t151 = t91 * t94;
t78 = t91 * qJD(3);
t13 = -t116 * t89 + t78;
t146 = qJD(4) * pkin(4);
t10 = t13 + t146;
t150 = -t13 + t10;
t149 = qJD(4) * t152 + t72 * t159;
t83 = t91 ^ 2;
t148 = t82 - t83;
t147 = t82 + t83;
t67 = t84 * pkin(2) + pkin(6);
t145 = qJ(5) + t67;
t46 = t84 * t90 - t86 * t92;
t143 = qJD(2) * t46;
t142 = qJD(2) * t91;
t141 = qJD(4) * t89;
t139 = qJDD(4) * pkin(4);
t47 = t84 * t92 + t86 * t90;
t38 = t47 * qJD(1);
t138 = t38 * qJD(2);
t137 = t143 * qJD(2);
t136 = t89 * qJD(3);
t29 = t86 * t60 - t63;
t70 = pkin(3) + t157;
t19 = -t70 * qJD(2) + qJD(5) - t29;
t134 = -qJD(5) - t19;
t133 = qJDD(1) - g(3);
t74 = t89 * qJDD(2);
t76 = t91 * qJDD(2);
t131 = qJ(5) * qJDD(2);
t129 = qJD(2) * qJD(4);
t128 = qJD(2) * qJD(5);
t126 = t168 * t84 - t86 * t41;
t124 = t147 * t40;
t122 = t89 * t129;
t121 = t91 * t129;
t5 = t89 * qJDD(3) - t22 * t141 + t169 * t91;
t120 = t91 * t138 + t40 * t141 + (g(1) * t153 + g(2) * t155) * t71;
t51 = -t70 - t158;
t8 = pkin(4) * t122 - t70 * qJDD(2) + qJDD(5) + t126;
t119 = qJDD(2) * t51 + t8;
t118 = qJD(4) * t145;
t117 = qJDD(2) * t147;
t115 = t89 * t121;
t14 = t136 + t111;
t113 = t10 * t89 - t14 * t91;
t17 = -t89 * t22 + t78;
t18 = t91 * t22 + t136;
t112 = t17 * t89 - t18 * t91;
t37 = t47 * qJD(2);
t110 = -t37 * qJD(2) - t46 * qJDD(2);
t11 = -qJDD(2) * pkin(3) + t126;
t68 = -pkin(3) - t158;
t93 = qJD(4) ^ 2;
t108 = qJDD(2) * t68 + t67 * t93 + t11;
t75 = t91 * qJDD(3);
t107 = -g(1) * (-t72 * t154 + t155) - g(2) * (-t72 * t156 - t153) + t71 * t159 + t75;
t105 = -t131 - t169;
t103 = t47 * t93 - t110;
t100 = 0.2e1 * t143 * qJD(4) - qJDD(4) * t47;
t21 = -qJD(2) * pkin(3) - t29;
t99 = -qJDD(4) * t67 + (qJD(2) * t68 + t21) * qJD(4);
t98 = -g(3) * t92 + t114 * t90;
t97 = -g(1) * (-t72 * t153 - t156) - g(2) * (-t72 * t155 + t154) - t5 + t91 * t161;
t96 = -t109 - t138;
t6 = -t18 * qJD(4) - t89 * t12 + t75;
t95 = t5 * t91 - t6 * t89 + (-t17 * t91 - t18 * t89) * qJD(4);
t88 = -qJ(5) - pkin(6);
t80 = t92 * pkin(2);
t64 = t89 * t151;
t54 = t148 * t94;
t53 = qJDD(4) * t91 - t93 * t89;
t52 = qJDD(4) * t89 + t93 * t91;
t45 = t83 * qJDD(2) - 0.2e1 * t115;
t44 = t82 * qJDD(2) + 0.2e1 * t115;
t43 = t145 * t91;
t42 = t145 * t89;
t25 = -t89 * qJD(5) - t91 * t118;
t24 = t91 * qJD(5) - t89 * t118;
t23 = -0.2e1 * t148 * t129 + 0.2e1 * t89 * t76;
t7 = t47 * t117 - t147 * t137;
t4 = t91 * t128 + (-t122 + t76) * qJ(5) + t5;
t3 = t139 + t75 - qJD(4) * t111 + (t105 - t128) * t89;
t2 = t100 * t89 - t103 * t91;
t1 = t100 * t91 + t103 * t89;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, 0, 0, 0, 0, 0, t92 * qJDD(2) - t94 * t90, -qJDD(2) * t90 - t94 * t92, 0, -g(3) + (t90 ^ 2 + t92 ^ 2) * qJDD(1), 0, 0, 0, 0, 0, 0, t110, -t47 * qJDD(2) + t137, 0, t126 * t46 - t143 * t30 + t16 * t47 - t29 * t37 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t7, t11 * t46 + t112 * t143 + t21 * t37 + t95 * t47 - g(3), 0, 0, 0, 0, 0, 0, t2, t1, t7, t19 * t37 + t8 * t46 - g(3) + t113 * t143 + (-t3 * t89 + t4 * t91 + (-t10 * t91 - t14 * t89) * qJD(4)) * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t77 + t98, t114 * t92 - t133 * t90, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t86 * t140 - t126 + t138 - t167, t40 * qJD(2) - t84 * t140 + t101 - t16, 0, t29 * t38 - t30 * t40 + (-t126 * t86 + t16 * t84 + t98) * pkin(2), t44, t23, t52, t45, t53, 0, t99 * t89 + (-t108 - t160) * t91 + t120, t99 * t91 + (t108 + t96) * t89 + t149, -qJD(2) * t124 + t67 * t117 - t101 + t95, t11 * t68 - t21 * t38 - g(3) * (t72 * pkin(3) + t71 * pkin(6) + t80) + t112 * t40 + t95 * t67 + t114 * (pkin(3) * t71 - pkin(6) * t72 + t166), t44, t23, t52, t45, t53, 0, -t42 * qJDD(4) + (-t119 - t160) * t91 + (t25 + (t19 + (t51 - t157) * qJD(2)) * t89) * qJD(4) + t120, -t43 * qJDD(4) + (t19 * t91 - t24 + (t51 * t91 + t165) * qJD(2)) * qJD(4) + (t119 + t96) * t89 + t149, (-qJD(4) * t10 + qJDD(2) * t43 + t4) * t91 + (-t14 * qJD(4) + qJDD(2) * t42 - t3) * t89 + (t24 * t91 - t25 * t89 - t124 + (t42 * t91 - t43 * t89) * qJD(4)) * qJD(2) - t101, t4 * t43 - t3 * t42 + t8 * t51 - g(3) * (t72 * t70 - t71 * t88 + t80) + (pkin(4) * t141 - t38) * t19 + (t24 - t152) * t14 + (t89 * t40 + t25) * t10 + t114 * (t70 * t71 + t72 * t88 + t166); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(3) + t125, 0, 0, 0, 0, 0, 0, t53, -t52, 0, -t112 * qJD(4) + t5 * t89 + t6 * t91 + t125, 0, 0, 0, 0, 0, 0, t53, -t52, 0, -t113 * qJD(4) + t3 * t91 + t4 * t89 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, t54, t74, t64, t76, qJDD(4), (-qJD(2) * t21 - t12) * t89 + t107, t17 * qJD(4) - t21 * t142 + t97, 0, 0, -t64, t54, t74, t64, t76, qJDD(4), 0.2e1 * t139 + (t14 - t111) * qJD(4) + (pkin(4) * t151 + t134 * qJD(2) + t105) * t89 + t107, -t94 * t165 - t91 * t131 + t13 * qJD(4) + (qJ(5) * t141 + t134 * t91) * qJD(2) + t97, -pkin(4) * t74 + (-t146 + t150) * t142, t150 * t14 + (t3 + t125 * t91 + (-t19 * qJD(2) + t101) * t89) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76 + 0.2e1 * t122, t74 + 0.2e1 * t121, -t147 * t94, t113 * qJD(2) + t167 + t8;];
tau_reg = t9;
