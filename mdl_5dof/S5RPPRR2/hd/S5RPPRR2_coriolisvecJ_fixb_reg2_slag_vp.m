% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:40:02
% DurationCPUTime: 1.43s
% Computational Cost: add. (2418->192), mult. (5253->253), div. (0->0), fcn. (3741->6), ass. (0->118)
t102 = sin(qJ(4));
t104 = cos(qJ(4));
t99 = cos(pkin(8));
t135 = t104 * t99;
t98 = sin(pkin(8));
t79 = -t102 * t98 + t135;
t101 = sin(qJ(5));
t103 = cos(qJ(5));
t128 = qJD(5) * t103;
t129 = qJD(5) * t101;
t78 = t102 * t99 + t104 * t98;
t134 = qJD(1) * t78;
t73 = t78 * qJD(4);
t60 = qJD(1) * t73;
t125 = qJD(1) * t135;
t133 = qJD(1) * t98;
t126 = t102 * t133;
t84 = qJD(4) * t126;
t61 = qJD(4) * t125 - t84;
t71 = t125 - t126;
t12 = t101 * t61 + t103 * t60 + t128 * t134 + t129 * t71;
t35 = t101 * t71 + t103 * t134;
t95 = qJD(4) + qJD(5);
t145 = t35 * t95;
t161 = -t12 + t145;
t113 = -t101 * t134 + t103 * t71;
t144 = t113 ^ 2;
t146 = t35 ^ 2;
t160 = t144 - t146;
t97 = qJD(1) * qJ(2);
t91 = qJD(3) + t97;
t82 = pkin(3) * t133 + t91;
t45 = pkin(4) * t134 + t82;
t159 = t45 * t35;
t143 = t35 * t113;
t158 = t79 * qJD(3);
t13 = qJD(5) * t113 - t101 * t60 + t103 * t61;
t147 = t113 * t95;
t157 = -t13 + t147;
t111 = -t101 * t78 + t103 * t79;
t130 = qJD(4) * t104;
t131 = qJD(4) * t102;
t74 = t130 * t99 - t131 * t98;
t106 = qJD(5) * t111 - t101 * t73 + t103 * t74;
t156 = t106 * t95;
t155 = t111 * t12;
t154 = t45 * t113;
t100 = -pkin(1) - qJ(3);
t153 = qJD(1) * t100;
t140 = t98 ^ 2 + t99 ^ 2;
t152 = qJD(3) * t140;
t110 = t158 * qJD(1);
t85 = qJD(2) + t153;
t122 = -pkin(6) * qJD(1) + t85;
t62 = t122 * t98;
t63 = t122 * t99;
t33 = t102 * t63 + t104 * t62;
t24 = -t33 * qJD(4) - t110;
t15 = t60 * pkin(7) + t24;
t27 = -pkin(7) * t134 + t33;
t120 = t101 * t15 - t129 * t27;
t109 = t78 * qJD(3);
t47 = t63 * t130;
t23 = -qJD(1) * t109 - t131 * t62 + t47;
t14 = -t61 * pkin(7) + t23;
t138 = t102 * t62;
t32 = t104 * t63 - t138;
t26 = -pkin(7) * t71 + t32;
t25 = qJD(4) * pkin(4) + t26;
t1 = (qJD(5) * t25 + t14) * t103 + t120;
t112 = t101 * t79 + t103 * t78;
t151 = t106 * t35 + t112 * t13;
t121 = -t101 * t14 + t103 * t15;
t136 = t103 * t27;
t6 = t101 * t25 + t136;
t2 = -qJD(5) * t6 + t121;
t150 = t1 * t112 + t106 * t6 + t111 * t2;
t149 = t71 ^ 2;
t96 = qJD(1) * qJD(2);
t127 = 0.2e1 * t96;
t148 = pkin(4) * t71;
t142 = t71 * t134;
t141 = -pkin(6) + t100;
t80 = t141 * t98;
t81 = t141 * t99;
t44 = t102 * t81 + t104 * t80;
t139 = t101 * t27;
t89 = pkin(3) * t98 + qJ(2);
t132 = t74 * qJD(4);
t123 = -pkin(4) * t95 - t25;
t46 = pkin(4) * t61 + t96;
t118 = -t101 * t74 - t103 * t73;
t43 = -t102 * t80 + t104 * t81;
t117 = qJD(1) * t140;
t115 = -t60 * t79 - t71 * t73;
t114 = t134 * t74 + t61 * t78;
t30 = -pkin(7) * t79 + t43;
t31 = -pkin(7) * t78 + t44;
t9 = -t101 * t31 + t103 * t30;
t10 = t101 * t30 + t103 * t31;
t108 = t23 * t78 + t24 * t79 - t32 * t73 + t33 * t74;
t28 = t130 * t81 - t131 * t80 - t109;
t29 = -qJD(4) * t44 - t158;
t105 = qJD(1) ^ 2;
t66 = t134 ^ 2;
t64 = t73 * qJD(4);
t54 = pkin(4) * t74 + qJD(2);
t52 = pkin(4) * t78 + t89;
t22 = t73 * pkin(7) + t29;
t21 = -t74 * pkin(7) + t28;
t18 = -qJD(5) * t112 + t118;
t17 = t128 * t78 + t129 * t79 - t118;
t8 = t103 * t26 - t139;
t7 = -t101 * t26 - t136;
t5 = t103 * t25 - t139;
t4 = -qJD(5) * t10 - t101 * t21 + t103 * t22;
t3 = qJD(5) * t9 + t101 * t22 + t103 * t21;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t127, qJ(2) * t127, 0, 0, 0, 0, 0, 0, t98 * t127, t99 * t127, 0.2e1 * qJD(3) * t117, (t91 + t97) * qJD(2) + (-t85 - t153) * t152, t115, t134 * t73 + t60 * t78 - t61 * t79 - t71 * t74, -t64, t114, -t132, 0, 0.2e1 * qJD(2) * t134 + t29 * qJD(4) + t89 * t61 + t82 * t74, -t28 * qJD(4) - t89 * t60 - t82 * t73 + (qJD(1) * t79 + t71) * qJD(2), -t134 * t28 - t29 * t71 + t43 * t60 - t44 * t61 - t108, t23 * t44 + t24 * t43 + t33 * t28 + t32 * t29 + (qJD(1) * t89 + t82) * qJD(2), -t113 * t17 - t155, -t106 * t113 - t111 * t13 + t112 * t12 + t17 * t35, -t17 * t95, t151, -t156, 0, t106 * t45 + t112 * t46 + t13 * t52 + t35 * t54 + t4 * t95, t111 * t46 + t113 * t54 - t12 * t52 - t17 * t45 - t3 * t95, -t10 * t13 - t113 * t4 + t12 * t9 + t17 * t5 - t3 * t35 - t150, t1 * t10 + t2 * t9 + t3 * t6 + t4 * t5 + t45 * t54 + t46 * t52; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t105 * qJ(2), 0, 0, 0, 0, 0, 0, -t105 * t98, -t105 * t99, 0, (-t91 - t152) * qJD(1), 0, 0, 0, 0, 0, 0, -qJD(1) * t134 - t64, -qJD(1) * t71 - t132, -t114 - t115, -qJD(1) * t82 + t108, 0, 0, 0, 0, 0, 0, -qJD(1) * t35 + t18 * t95, -qJD(1) * t113 - t156, -t113 * t18 - t151 + t155, -qJD(1) * t45 + t18 * t5 + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140 * t105, t117 * t85 + t96, 0, 0, 0, 0, 0, 0, -t84 + (t71 + t125) * qJD(4), -0.2e1 * t134 * qJD(4), -t66 - t149, t134 * t33 + t32 * t71 + t96, 0, 0, 0, 0, 0, 0, t13 + t147, -t12 - t145, -t144 - t146, t113 * t5 + t35 * t6 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t142, -t66 + t149, 0, -t142, t84 + (t71 - t125) * qJD(4), 0, -t82 * t71 - t110, -t47 + (t32 + t138) * qJD(4) + (t82 + qJD(3)) * t134, 0, 0, t143, t160, t161, -t143, t157, 0, -t35 * t148 - t154 - t7 * t95 + (t101 * t123 - t136) * qJD(5) + t121, -t113 * t148 + t159 + t8 * t95 + (qJD(5) * t123 - t14) * t103 - t120, t6 * t113 + t8 * t35 - t5 * t35 + t7 * t113 + (-t101 * t13 + t103 * t12 + (t101 * t113 - t103 * t35) * qJD(5)) * pkin(4), -t5 * t7 - t6 * t8 + (t1 * t101 + t103 * t2 - t45 * t71 + (-t101 * t5 + t103 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t143, t160, t161, -t143, t157, 0, t6 * t95 - t154 + t2, t5 * t95 - t1 + t159, 0, 0;];
tauc_reg = t11;
