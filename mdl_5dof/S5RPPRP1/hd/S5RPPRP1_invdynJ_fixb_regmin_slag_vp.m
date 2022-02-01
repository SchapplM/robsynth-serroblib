% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S5RPPRP1
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tau_reg [5x18]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S5RPPRP1_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_invdynJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:49
% DurationCPUTime: 1.21s
% Computational Cost: add. (1100->217), mult. (2193->294), div. (0->0), fcn. (1410->10), ass. (0->133)
t68 = qJ(1) + pkin(7);
t63 = sin(t68);
t64 = cos(t68);
t106 = -g(1) * t63 + g(2) * t64;
t74 = cos(pkin(7));
t59 = -t74 * pkin(1) - pkin(2);
t123 = qJDD(1) * t59;
t42 = qJDD(3) + t123;
t164 = -t106 - t42;
t71 = sin(pkin(8));
t117 = qJD(1) * qJD(4);
t78 = cos(qJ(4));
t105 = t78 * t117;
t76 = sin(qJ(4));
t85 = qJDD(1) * t76 + t105;
t163 = t85 * t71;
t73 = cos(pkin(8));
t126 = t73 * qJD(1);
t51 = -qJD(4) + t126;
t162 = qJD(4) + t51;
t118 = qJD(1) * qJD(3);
t72 = sin(pkin(7));
t56 = t72 * pkin(1) + qJ(3);
t39 = qJDD(1) * t56 + t118;
t161 = -t71 * (-qJ(5) - pkin(6)) + (t78 * pkin(4) + pkin(3)) * t73;
t143 = t73 * t76;
t26 = t63 * t143 + t64 * t78;
t28 = -t64 * t143 + t63 * t78;
t160 = -g(1) * t28 + g(2) * t26;
t152 = g(3) * t76;
t159 = t71 * t152 + t160;
t36 = -t73 * pkin(3) - t71 * pkin(6) + t59;
t25 = t36 * qJD(1) + qJD(3);
t44 = t56 * qJD(1);
t32 = t71 * qJD(2) + t73 * t44;
t104 = t78 * t25 - t76 * t32;
t135 = qJD(1) * t71;
t108 = qJ(5) * t135;
t6 = -t78 * t108 + t104;
t3 = -t51 * pkin(4) + t6;
t158 = -t6 + t3;
t156 = pkin(4) * t76;
t120 = t73 * qJDD(1);
t50 = -qJDD(4) + t120;
t151 = t50 * pkin(4);
t145 = t71 * t39;
t61 = t73 * qJDD(2);
t23 = -t61 + t145;
t150 = t23 * t71;
t149 = t50 * t73;
t147 = t64 * t76;
t66 = t71 ^ 2;
t80 = qJD(1) ^ 2;
t146 = t66 * t80;
t142 = t73 * t78;
t128 = qJD(4) * t78;
t131 = qJD(3) * t78;
t141 = t36 * t128 + t73 * t131;
t40 = t56 * t142;
t140 = t76 * t36 + t40;
t139 = t73 ^ 2 + t66;
t69 = t76 ^ 2;
t70 = t78 ^ 2;
t138 = -t69 - t70;
t137 = t69 - t70;
t136 = qJ(5) * t71;
t134 = qJD(1) * t76;
t133 = qJD(1) * t78;
t132 = qJD(3) * t73;
t130 = qJD(4) * t32;
t129 = qJD(4) * t76;
t127 = qJD(5) * t71;
t62 = t73 * qJD(2);
t15 = qJD(5) - t62 + (pkin(4) * t134 + t44) * t71;
t125 = qJD(5) + t15;
t121 = qJDD(1) * t78;
t119 = qJ(5) * qJDD(1);
t116 = qJD(1) * qJD(5);
t115 = t76 * t146;
t79 = cos(qJ(1));
t114 = t79 * pkin(1) + t64 * pkin(2) + t63 * qJ(3);
t113 = t78 * t136;
t112 = t71 * t134;
t111 = t56 * t129;
t110 = t51 * t129;
t31 = t71 * t44 - t62;
t109 = t31 * t135;
t77 = sin(qJ(1));
t107 = -t77 * pkin(1) + t64 * qJ(3);
t22 = t36 * qJDD(1) + qJDD(3);
t24 = t71 * qJDD(2) + t73 * t39;
t103 = t25 * t128 - t32 * t129 + t76 * t22 + t78 * t24;
t35 = (t56 + t156) * t71;
t102 = qJD(1) * t35 + t15;
t101 = -qJD(4) * t25 - t24;
t100 = t50 - t120;
t99 = t50 + t120;
t98 = pkin(4) * t163 + qJDD(5) - t61;
t97 = qJD(4) * t112;
t96 = -g(1) * t26 - g(2) * t28;
t27 = -t63 * t142 + t147;
t29 = t64 * t142 + t63 * t76;
t95 = -g(1) * t27 - g(2) * t29;
t94 = -g(1) * t64 - g(2) * t63;
t92 = g(1) * t77 - g(2) * t79;
t88 = -t76 * t25 - t78 * t32;
t7 = -t76 * t108 - t88;
t91 = t3 * t78 + t7 * t76;
t90 = t3 * t76 - t7 * t78;
t89 = t24 * t73 + t150;
t87 = t31 * t71 + t32 * t73;
t13 = t98 + t145;
t38 = (pkin(4) * t128 + qJD(3)) * t71;
t86 = qJD(1) * t38 + qJDD(1) * t35 + t13;
t18 = t78 * t22;
t84 = qJ(5) * t97 + t101 * t76 + t18;
t83 = -t51 ^ 2 - t146;
t82 = g(3) * t71 * t78 + g(1) * t29 - g(2) * t27 - t103;
t48 = t71 * t121;
t41 = t73 * t97;
t37 = t71 * t51 * t133;
t34 = t78 * t36;
t14 = -t76 * t136 + t140;
t12 = -t113 + t34 + (-t56 * t76 - pkin(4)) * t73;
t11 = t76 * t50 + t83 * t78;
t10 = -t78 * t50 + t83 * t76;
t9 = (t100 * t76 + (t51 - t126) * t128) * t71;
t8 = t41 + (t100 * t78 - t110) * t71;
t5 = -t76 * t132 - t78 * t127 + (-t40 + (-t36 + t136) * t76) * qJD(4);
t4 = -t76 * t127 + (-t56 * t143 - t113) * qJD(4) + t141;
t2 = (-t85 * qJ(5) - t76 * t116) * t71 + t103;
t1 = -t151 + (-t130 + (-t116 - t119) * t71) * t78 + t84;
t16 = [qJDD(1), t92, g(1) * t79 + g(2) * t77, (t92 + (t72 ^ 2 + t74 ^ 2) * qJDD(1) * pkin(1)) * pkin(1), (-t123 + t164) * t73, t39 * t139 + t89 + t94, t42 * t59 - g(1) * (-t63 * pkin(2) + t107) - g(2) * t114 + t89 * t56 + t87 * qJD(3), (qJDD(1) * t70 - 0.2e1 * t76 * t105) * t66, 0.2e1 * (t137 * t117 - t76 * t121) * t66, t41 + (-t99 * t78 + t110) * t71, (t99 * t76 + (t51 + t126) * t128) * t71, t149, -t18 * t73 - t34 * t50 + ((qJD(1) * t66 + t51 * t73) * t56 + t87) * t128 + (-(-qJD(4) * t36 - t132) * t51 - t101 * t73 + t66 * t118 + t150 + (t66 * qJDD(1) + t149) * t56) * t76 + t95, (-t73 * t111 + t141) * t51 + t140 * t50 + t103 * t73 + (-t31 * t129 + t23 * t78) * t71 + (t56 * t121 + (-t111 + t131) * qJD(1)) * t66 + t96, -t1 * t73 - t12 * t50 - t5 * t51 + (t102 * t128 + t86 * t76) * t71 + t95, t14 * t50 + t2 * t73 + t4 * t51 + (-t102 * t129 + t86 * t78) * t71 + t96, ((-qJD(4) * t7 - qJDD(1) * t12 - t1 + (-qJD(4) * t14 - t5) * qJD(1)) * t78 + (qJD(4) * t3 - qJDD(1) * t14 - t2 + (qJD(4) * t12 - t4) * qJD(1)) * t76 - t106) * t71, t2 * t14 + t7 * t4 + t1 * t12 + t3 * t5 + t13 * t35 + t15 * t38 - g(1) * (pkin(4) * t147 + t107) - g(2) * (t161 * t64 + t114) + (-g(1) * (-pkin(2) - t161) - g(2) * t156) * t63; 0, 0, 0, qJDD(2) - g(3), 0, 0, -t23 * t73 + t24 * t71 - g(3), 0, 0, 0, 0, 0, t9, t8, t9, t8, 0, -t13 * t73 - g(3) + (-t91 * qJD(4) - t1 * t76 + t2 * t78) * t71; 0, 0, 0, 0, -t120, -t139 * t80, -t87 * qJD(1) - t164, 0, 0, 0, 0, 0, t10, t11, t10, t11, t138 * t71 * qJDD(1), t1 * t78 + t2 * t76 - t90 * qJD(4) + (-t15 * t71 + t90 * t73) * qJD(1) + t106; 0, 0, 0, 0, 0, 0, 0, t78 * t115, -t137 * t146, -t162 * t112 + t48, -t37 - t163, -t50, -t78 * t109 + t162 * t88 - t76 * t24 + t159 + t18, -t104 * t51 + t76 * t109 + t82, -0.2e1 * t151 - t7 * t51 + (-pkin(4) * t115 - t130 + (-t125 * qJD(1) - t119) * t71) * t78 + t84 + t159, -t70 * pkin(4) * t146 - t6 * t51 + (t76 * t119 + (qJ(5) * t128 + t125 * t76) * qJD(1)) * t71 + t82, (-pkin(4) * t121 + (pkin(4) * qJD(4) - t158) * t134) * t71, t158 * t7 + (t1 + (-t15 * t133 + t152) * t71 + t160) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37 + t163, t48 + (-qJD(4) + t51) * t112, t138 * t146, g(3) * t73 + (t91 * qJD(1) + t39 + t94) * t71 + t98;];
tau_reg = t16;
