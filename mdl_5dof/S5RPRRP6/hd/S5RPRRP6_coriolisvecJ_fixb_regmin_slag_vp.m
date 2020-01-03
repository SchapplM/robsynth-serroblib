% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% tauc_reg [5x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:43:16
% EndTime: 2019-12-31 18:43:19
% DurationCPUTime: 0.90s
% Computational Cost: add. (1093->193), mult. (2592->281), div. (0->0), fcn. (1538->6), ass. (0->115)
t63 = sin(pkin(8)) * pkin(1) + pkin(6);
t57 = t63 * qJD(1);
t74 = sin(qJ(3));
t76 = cos(qJ(3));
t152 = t76 * qJD(2) - t74 * t57;
t73 = sin(qJ(4));
t117 = qJD(4) * t73;
t105 = t74 * t117;
t75 = cos(qJ(4));
t113 = t75 * qJD(3);
t111 = t76 * qJD(1);
t62 = -qJD(4) + t111;
t69 = t74 ^ 2;
t86 = qJD(1) * t69 - t62 * t76;
t151 = -t62 * t105 - t86 * t113;
t116 = qJD(4) * t75;
t104 = t74 * t116;
t108 = qJD(3) * qJD(4);
t114 = t73 * qJD(3);
t20 = (t76 * t114 + t104) * qJD(1) + t73 * t108;
t121 = qJD(1) * t74;
t52 = t75 * t121 + t114;
t150 = t52 ^ 2;
t67 = t74 * qJD(2);
t31 = t76 * t57 + t67;
t24 = qJD(3) * pkin(7) + t31;
t64 = -cos(pkin(8)) * pkin(1) - pkin(2);
t44 = -t76 * pkin(3) - t74 * pkin(7) + t64;
t27 = t44 * qJD(1);
t9 = -t73 * t24 + t75 * t27;
t6 = -t52 * qJ(5) + t9;
t5 = -t62 * pkin(4) + t6;
t149 = t5 - t6;
t122 = t75 * qJ(5);
t83 = pkin(4) * t74 - t76 * t122;
t131 = -qJ(5) - pkin(7);
t93 = qJD(4) * t131;
t89 = pkin(3) * t74 - pkin(7) * t76;
t55 = t89 * qJD(1);
t94 = -t152 * t73 + t75 * t55;
t148 = -t83 * qJD(1) - t73 * qJD(5) + t75 * t93 - t94;
t23 = -qJD(3) * pkin(3) - t152;
t147 = t23 * t73;
t146 = t23 * t75;
t119 = qJD(3) * t76;
t26 = qJD(3) * t67 + t57 * t119;
t145 = t26 * t73;
t144 = t26 * t75;
t102 = t73 * t121;
t50 = t102 - t113;
t143 = t50 * t62;
t142 = t50 * t74;
t141 = t52 * t62;
t140 = t62 * t73;
t139 = t62 * t75;
t138 = t73 * t27;
t137 = t73 * t76;
t136 = t74 * t75;
t135 = t75 * t76;
t134 = t76 * t20;
t77 = qJD(3) ^ 2;
t133 = t77 * t74;
t132 = t77 * t76;
t101 = t76 * t113;
t130 = -t50 * t101 - t20 * t136;
t129 = t152 * t75 + t73 * t55;
t56 = t89 * qJD(3);
t128 = t44 * t116 + t73 * t56;
t112 = t75 * qJD(5);
t123 = t73 * qJ(5);
t127 = t111 * t123 + t73 * t93 + t112 - t129;
t126 = t74 * t63 * t114 + t75 * t56;
t54 = t63 * t135;
t125 = t73 * t44 + t54;
t124 = -t76 ^ 2 + t69;
t58 = qJD(1) * t64;
t120 = qJD(3) * t74;
t118 = qJD(4) * t50;
t115 = t23 * qJD(4);
t109 = qJD(1) * qJD(3);
t25 = t152 * qJD(3);
t43 = qJD(1) * t56;
t107 = t27 * t116 + t75 * t25 + t73 * t43;
t106 = t52 * t119;
t103 = t62 * t116;
t100 = pkin(4) * t73 + t63;
t98 = t74 * t109;
t19 = -qJD(1) * t101 + qJD(4) * t102 - t75 * t108;
t97 = t52 * t120 + t19 * t76;
t96 = t62 * t63 + t24;
t95 = t73 * t25 - t75 * t43;
t91 = t74 * t103;
t90 = t52 * t104;
t11 = t20 * pkin(4) + t26;
t10 = t75 * t24 + t138;
t7 = -t50 * qJ(5) + t10;
t88 = -t5 * t75 - t7 * t73;
t87 = t5 * t73 - t7 * t75;
t84 = 0.2e1 * qJD(3) * t58;
t82 = -t24 * t117 + t107;
t81 = t86 * t73;
t79 = -t10 * qJD(4) - t95;
t78 = qJD(1) ^ 2;
t60 = t131 * t75;
t59 = t131 * t73;
t47 = t50 ^ 2;
t35 = t75 * t44;
t15 = t50 * pkin(4) + qJD(5) + t23;
t14 = -t74 * t123 + t125;
t13 = -t74 * t122 + t35 + (-t63 * t73 - pkin(4)) * t76;
t4 = (-qJ(5) * qJD(4) - qJD(3) * t63) * t136 + (-qJD(5) * t74 + (-qJ(5) * qJD(3) - qJD(4) * t63) * t76) * t73 + t128;
t3 = -t74 * t112 + t83 * qJD(3) + (-t54 + (qJ(5) * t74 - t44) * t73) * qJD(4) + t126;
t2 = -t20 * qJ(5) - t50 * qJD(5) + t82;
t1 = pkin(4) * t98 + t19 * qJ(5) - t52 * qJD(5) + t79;
t8 = [0, 0, 0, 0, 0.2e1 * t76 * t98, -0.2e1 * t124 * t109, t132, -t133, 0, -t63 * t132 + t74 * t84, t63 * t133 + t76 * t84, -t19 * t136 + (t101 - t105) * t52, -t90 + (-t106 + (t19 + t118) * t74) * t73 + t130, -t151 + t97, t91 + t134 + (-t81 - t142) * qJD(3), (-t62 - t111) * t120, -(-t44 * t117 + t126) * t62 + ((t50 * t63 + t147) * qJD(3) + (t96 * t75 + t138) * qJD(4) + t95) * t76 + (t75 * t115 + t63 * t20 + t145 + ((-t63 * t137 + t35) * qJD(1) + t9) * qJD(3)) * t74, t128 * t62 + (-t96 * t117 + (t52 * t63 + t146) * qJD(3) + t107) * t76 + (-t73 * t115 - t63 * t19 + t144 + (-t125 * qJD(1) - t63 * t139 - t10) * qJD(3)) * t74, t13 * t19 - t14 * t20 - t3 * t52 - t4 * t50 + t88 * t119 + (t87 * qJD(4) - t1 * t75 - t2 * t73) * t74, t1 * t13 + t2 * t14 + t5 * t3 + t7 * t4 + t15 * t100 * t119 + (t15 * pkin(4) * t116 + t11 * t100) * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t133, -t132, 0, 0, 0, 0, 0, t91 - t134 + (-t81 + t142) * qJD(3), t151 + t97, t90 + (t106 + (-t19 + t118) * t74) * t73 + t130, (-t87 * qJD(3) - t11) * t76 + (qJD(3) * t15 + t88 * qJD(4) - t1 * t73 + t2 * t75) * t74; 0, 0, 0, 0, -t74 * t78 * t76, t124 * t78, 0, 0, 0, t31 * qJD(3) - t58 * t121 - t26, -t58 * t111, -t52 * t139 - t19 * t73, (-t19 + t143) * t75 + (-t20 + t141) * t73, -t103 + (t62 * t135 + (-t52 + t114) * t74) * qJD(1), t62 * t117 + (-t62 * t137 + (t50 + t113) * t74) * qJD(1), t62 * t121, -pkin(3) * t20 - t144 + t94 * t62 - t31 * t50 + (pkin(7) * t139 + t147) * qJD(4) + (-t9 * t74 + (-pkin(7) * t120 - t23 * t76) * t73) * qJD(1), pkin(3) * t19 + t145 - t129 * t62 - t31 * t52 + (-pkin(7) * t140 + t146) * qJD(4) + (-t23 * t135 + (-pkin(7) * t113 + t10) * t74) * qJD(1), t59 * t19 + t60 * t20 - t148 * t52 - t127 * t50 + (t62 * t5 + t2) * t75 + (t62 * t7 - t1) * t73, -t2 * t60 + t1 * t59 + t11 * (-t75 * pkin(4) - pkin(3)) + t127 * t7 + t148 * t5 + (-pkin(4) * t140 - t31) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t50, -t47 + t150, -t19 - t143, -t141 - t20, t98, -t10 * t62 - t23 * t52 + t79, t23 * t50 - t9 * t62 - t82, pkin(4) * t19 - t149 * t50, t149 * t7 + (-t15 * t52 + t1) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47 - t150, t5 * t52 + t7 * t50 + t11;];
tauc_reg = t8;
