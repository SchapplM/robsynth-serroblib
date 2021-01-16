% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauc_reg [5x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:53
% EndTime: 2021-01-15 19:26:01
% DurationCPUTime: 1.43s
% Computational Cost: add. (1439->236), mult. (3052->344), div. (0->0), fcn. (1717->4), ass. (0->130)
t67 = cos(qJ(4));
t107 = t67 * qJD(3);
t68 = cos(qJ(3));
t115 = qJD(1) * t68;
t65 = sin(qJ(4));
t45 = t65 * t115 - t107;
t66 = sin(qJ(3));
t108 = t66 * qJD(1);
t59 = qJD(4) + t108;
t142 = t45 * t59;
t94 = t66 * t107;
t110 = qJD(4) * t68;
t97 = t65 * t110;
t75 = t94 + t97;
t18 = t75 * qJD(1) - qJD(4) * t107;
t154 = -t18 - t142;
t100 = 0.2e1 * qJD(1);
t69 = -pkin(1) - pkin(6);
t132 = t66 * t69;
t52 = t66 * pkin(3) - t68 * pkin(7) + qJ(2);
t121 = t67 * t132 + t65 * t52;
t109 = t65 * qJD(3);
t93 = t67 * t115;
t47 = t93 + t109;
t141 = t47 * t59;
t95 = t66 * t109;
t19 = -qJD(1) * t95 + t47 * qJD(4);
t153 = t19 + t141;
t152 = t47 ^ 2;
t33 = t52 * qJD(1);
t58 = t69 * qJD(1) + qJD(2);
t51 = t66 * t58;
t36 = qJD(3) * pkin(7) + t51;
t13 = t67 * t33 - t65 * t36;
t8 = -t47 * qJ(5) + t13;
t7 = t59 * pkin(4) + t8;
t151 = t7 - t8;
t150 = t45 * pkin(4);
t14 = t65 * t33 + t67 * t36;
t9 = -t45 * qJ(5) + t14;
t149 = t9 * t59;
t114 = qJD(3) * t66;
t11 = t19 * pkin(4) + t58 * t114;
t148 = t11 * t65;
t147 = t11 * t67;
t146 = t18 * t65;
t145 = t18 * t66;
t144 = t19 * t66;
t128 = t68 * t58;
t37 = -qJD(3) * pkin(3) - t128;
t143 = t37 * t65;
t140 = t47 * t65;
t139 = t58 * t65;
t138 = t59 * t67;
t137 = t65 * t59;
t136 = t65 * t66;
t135 = t65 * t68;
t134 = t65 * t69;
t133 = t66 * t67;
t131 = t67 * t68;
t130 = t68 * t18;
t129 = t68 * t19;
t70 = qJD(3) ^ 2;
t127 = t70 * t66;
t126 = t70 * t68;
t125 = -qJ(5) - pkin(7);
t102 = t65 * t128;
t117 = t67 * qJ(5);
t82 = pkin(3) * t68 + pkin(7) * t66;
t50 = t82 * qJD(1);
t35 = t67 * t50;
t88 = qJD(4) * t125;
t124 = -t102 + t35 + (pkin(4) * t68 + t66 * t117) * qJD(1) + t65 * qJD(5) - t67 * t88;
t106 = t67 * qJD(5);
t122 = t67 * t128 + t65 * t50;
t99 = t65 * t108;
t123 = qJ(5) * t99 - t65 * t88 - t106 + t122;
t64 = t68 ^ 2;
t120 = t66 ^ 2 - t64;
t71 = qJD(1) ^ 2;
t119 = -t70 - t71;
t118 = pkin(4) * qJD(1);
t116 = t71 * qJ(2);
t113 = qJD(3) * t68;
t112 = qJD(4) * t65;
t111 = qJD(4) * t67;
t105 = qJ(2) * qJD(3);
t104 = qJD(1) * qJD(3);
t103 = t65 * t132;
t43 = t82 * qJD(3) + qJD(2);
t98 = t68 * t107;
t101 = t52 * t111 + t65 * t43 + t69 * t98;
t96 = t67 * t110;
t92 = qJD(2) * t100;
t91 = pkin(4) - t134;
t90 = t68 * t104;
t89 = -qJD(5) - t150;
t27 = t43 * qJD(1);
t87 = -t33 * t111 + t36 * t112 - t65 * t27 - t58 * t98;
t86 = t59 + t108;
t85 = t45 + t107;
t84 = -t47 + t109;
t83 = qJD(4) * t66 + qJD(1);
t81 = t65 * t9 + t67 * t7;
t80 = t65 * t7 - t67 * t9;
t79 = qJD(1) * t64 - t59 * t66;
t77 = -pkin(7) * t113 + t37 * t66;
t76 = t19 * qJ(5) + t87;
t22 = t67 * t27;
t74 = -t14 * qJD(4) + t22;
t72 = t18 * qJ(5) + t74;
t61 = -t67 * pkin(4) - pkin(3);
t55 = t125 * t67;
t54 = t125 * t65;
t44 = (pkin(4) * t65 - t69) * t68;
t42 = t45 ^ 2;
t41 = t67 * t52;
t31 = t67 * t43;
t23 = -pkin(4) * t99 + t51;
t20 = t69 * t114 + (-t95 + t96) * pkin(4);
t17 = -qJ(5) * t135 + t121;
t16 = t37 - t89;
t15 = -t68 * t117 + t91 * t66 + t41;
t6 = -qJ(5) * t96 + (-qJD(5) * t68 + (qJ(5) * qJD(3) - qJD(4) * t69) * t66) * t65 + t101;
t5 = qJ(5) * t94 + t31 - t121 * qJD(4) + (qJ(5) * t112 + t91 * qJD(3) - t106) * t68;
t4 = -t129 - t83 * t138 + (-t86 * t135 + t66 * t45) * qJD(3);
t3 = t130 + t83 * t137 + (-t59 * t131 + (t47 - t93) * t66) * qJD(3);
t2 = -t45 * qJD(5) - t76;
t1 = -t47 * qJD(5) + (t118 - t139) * t113 + t72;
t10 = [0, 0, 0, 0, t92, qJ(2) * t92, -0.2e1 * t66 * t90, 0.2e1 * t120 * t104, -t127, -t126, 0, -t69 * t127 + (qJD(2) * t66 + t68 * t105) * t100, -t69 * t126 + (qJD(2) * t68 - t66 * t105) * t100, -t67 * t130 - t47 * t75, (t45 * t67 + t140) * t114 + (t146 - t19 * t67 + (t45 * t65 - t47 * t67) * qJD(4)) * t68, -t59 * t97 - t145 + (t47 * t68 + t67 * t79) * qJD(3), -t59 * t96 - t144 + (-t45 * t68 - t65 * t79) * qJD(3), t86 * t113, -t69 * t129 + t22 * t66 + t31 * t59 + (-t121 * t59 + t37 * t131 - t14 * t66) * qJD(4) + ((t69 * t45 - t143) * t66 + (-t59 * t134 + (t41 - t103) * qJD(1) + t13) * t68) * qJD(3), -(-qJD(4) * t103 + t101) * t59 + t87 * t66 + (-t37 * t112 + t69 * t18) * t68 + ((-t121 * qJD(1) - t14) * t68 + (t69 * t47 + (-t37 + t128) * t67) * t66) * qJD(3), t44 * t19 + t20 * t45 + t5 * t59 + (-t16 * t109 + t1) * t66 + (t16 * t111 + t148 + (qJD(1) * t15 + t7) * qJD(3)) * t68, -t44 * t18 + t20 * t47 - t6 * t59 + (-t16 * t107 - t2) * t66 + (-t16 * t112 + t147 + (-qJD(1) * t17 - t9) * qJD(3)) * t68, t15 * t18 - t17 * t19 - t6 * t45 - t5 * t47 + t81 * t114 + (qJD(4) * t80 - t1 * t67 - t2 * t65) * t68, t1 * t15 + t11 * t44 + t16 * t20 + t2 * t17 + t7 * t5 + t9 * t6; 0, 0, 0, 0, -t71, -t116, 0, 0, 0, 0, 0, t119 * t66, t119 * t68, 0, 0, 0, 0, 0, t4, t3, t4, t3, (-t45 * t113 + t47 * t83 - t144) * t67 + (t47 * t113 + t45 * t83 - t145) * t65, -t81 * qJD(1) + (-qJD(3) * t80 - t11) * t68 + (qJD(3) * t16 - qJD(4) * t81 - t1 * t65 + t2 * t67) * t66; 0, 0, 0, 0, 0, 0, t68 * t71 * t66, -t120 * t71, 0, 0, 0, -t68 * t116, t66 * t116, t47 * t138 - t146, -t153 * t65 + t154 * t67, t59 * t111 + (t59 * t133 + t68 * t84) * qJD(1), -t59 * t112 + (-t59 * t136 + t68 * t85) * qJD(1), -t59 * t115, -pkin(3) * t19 - t35 * t59 + (t59 * t135 - t66 * t85) * t58 + (-pkin(7) * t138 + t143) * qJD(4) + (-t13 * t68 + t65 * t77) * qJD(1), pkin(3) * t18 + t122 * t59 + t84 * t51 + (pkin(7) * t137 + t37 * t67) * qJD(4) + (t14 * t68 + t67 * t77) * qJD(1), -t147 + t61 * t19 - t23 * t45 - t124 * t59 + (t16 + t150) * t112 + (t16 * t136 + (qJD(3) * t54 - t7) * t68) * qJD(1), t148 - t61 * t18 - t23 * t47 + t123 * t59 + (pkin(4) * t140 + t16 * t67) * qJD(4) + (t16 * t133 + (qJD(3) * t55 + t9) * t68) * qJD(1), t54 * t18 + t55 * t19 + t124 * t47 + t123 * t45 + (-t59 * t7 + t2) * t67 + (-t1 - t149) * t65, t1 * t54 + t11 * t61 - t2 * t55 - t123 * t9 - t124 * t7 + (pkin(4) * t112 - t23) * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47 * t45, -t42 + t152, -t18 + t142, t141 - t19, t90, -qJD(3) * t102 + t14 * t59 - t37 * t47 + t74, t13 * t59 + t37 * t45 + t87, t149 + (0.2e1 * t118 - t139) * t113 + (-t16 + t89) * t47 + t72, -t152 * pkin(4) + t8 * t59 + (qJD(5) + t16) * t45 + t76, t18 * pkin(4) - t151 * t45, t151 * t9 + (-t16 * t47 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t153, t154, -t42 - t152, t9 * t45 + t7 * t47 + t11;];
tauc_reg = t10;
