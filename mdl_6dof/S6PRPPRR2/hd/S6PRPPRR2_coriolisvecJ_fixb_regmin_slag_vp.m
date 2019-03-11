% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
% 
% Output:
% tauc_reg [6x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PRPPRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:20:01
% EndTime: 2019-03-08 19:20:04
% DurationCPUTime: 0.83s
% Computational Cost: add. (828->187), mult. (2089->285), div. (0->0), fcn. (1566->10), ass. (0->122)
t56 = sin(pkin(6));
t110 = qJD(1) * t56;
t64 = cos(qJ(2));
t91 = t64 * t110;
t43 = qJD(2) * pkin(2) + t91;
t55 = sin(pkin(11));
t57 = cos(pkin(11));
t61 = sin(qJ(2));
t95 = t61 * t110;
t22 = t55 * t43 + t57 * t95;
t16 = qJD(2) * qJ(4) + t22;
t32 = (t55 * t64 + t57 * t61) * t56;
t27 = qJD(1) * t32;
t135 = -t16 + t27;
t60 = sin(qJ(5));
t99 = t60 * qJD(2);
t49 = qJD(6) + t99;
t134 = -qJD(6) + t49;
t63 = cos(qJ(5));
t54 = t63 ^ 2;
t109 = qJD(2) * t54;
t59 = sin(qJ(6));
t103 = qJD(6) * t59;
t93 = t63 * t103;
t62 = cos(qJ(6));
t98 = t62 * qJD(5);
t133 = -(-t49 * t60 + t109) * t98 + t49 * t93;
t100 = t59 * qJD(5);
t107 = qJD(2) * t63;
t89 = t62 * t107;
t39 = t89 + t100;
t90 = t60 * t100;
t26 = -qJD(2) * t90 + t39 * qJD(6);
t44 = t55 * t95;
t21 = t57 * t43 - t44;
t77 = qJD(4) - t21;
t14 = (-pkin(3) - pkin(8)) * qJD(2) + t77;
t58 = cos(pkin(6));
t47 = t58 * qJD(1) + qJD(3);
t10 = t60 * t14 + t63 * t47;
t28 = qJD(2) * t32;
t23 = qJD(1) * t28;
t4 = t10 * qJD(5) - t63 * t23;
t132 = t4 * t59;
t131 = t4 * t62;
t130 = t55 * pkin(2);
t122 = t57 * t64;
t124 = t55 * t61;
t31 = (-t122 + t124) * t56;
t129 = t23 * t31;
t69 = -t60 * t98 - t93;
t25 = t69 * qJD(2) + qJD(6) * t98;
t128 = t25 * t59;
t37 = t59 * t107 - t98;
t127 = t37 * t49;
t126 = t39 * t49;
t125 = t49 * t59;
t66 = qJD(2) ^ 2;
t123 = t56 * t66;
t121 = t59 * t60;
t120 = t60 * t26;
t119 = t62 * t49;
t118 = t63 * t25;
t117 = t63 * t37;
t65 = qJD(5) ^ 2;
t116 = t65 * t60;
t115 = t65 * t63;
t104 = qJD(5) * t63;
t114 = t39 * t104 + t25 * t60;
t81 = t57 * t91;
t30 = -t44 + t81;
t79 = pkin(5) * t63 + pkin(9) * t60;
t36 = t79 * qJD(5) + qJD(4);
t113 = t30 - t36;
t108 = qJD(2) * t56;
t82 = t108 * t124;
t24 = -qJD(1) * t82 + qJD(2) * t81;
t112 = t60 ^ 2 - t54;
t111 = -t65 - t66;
t88 = -t57 * pkin(2) - pkin(3);
t48 = -pkin(8) + t88;
t106 = qJD(5) * t48;
t105 = qJD(5) * t60;
t102 = qJD(6) * t62;
t101 = t32 * qJD(5);
t97 = qJD(4) - t30;
t96 = qJD(2) * qJD(5);
t94 = t59 * t109;
t92 = t49 * t102;
t8 = qJD(5) * pkin(9) + t10;
t87 = t48 * t49 + t8;
t86 = t63 * t96;
t85 = t16 * qJD(2) - t23;
t84 = t49 + t99;
t80 = qJD(6) * t60 + qJD(2);
t50 = qJ(4) + t130;
t78 = qJD(2) * t50 - t135;
t71 = t60 * pkin(5) - t63 * pkin(9) + qJ(4);
t13 = t71 * qJD(2) + t22;
t1 = t62 * t13 - t59 * t8;
t2 = t59 * t13 + t62 * t8;
t9 = t63 * t14 - t60 * t47;
t18 = t31 * t60 + t58 * t63;
t76 = t18 * t62 + t32 * t59;
t75 = -t18 * t59 + t32 * t62;
t74 = t31 * t63 - t58 * t60;
t72 = t49 * t90 - t63 * t92;
t7 = -qJD(5) * pkin(5) - t9;
t70 = -pkin(9) * t104 + t60 * t7;
t19 = qJD(4) * qJD(2) + t24;
t68 = t97 * qJD(2) - t48 * t65 + t19;
t3 = t9 * qJD(5) + t60 * t23;
t67 = -qJD(5) * t7 - qJD(6) * t13 + t27 * t49 - t3;
t42 = t79 * qJD(2);
t34 = t71 + t130;
t29 = t108 * t122 - t82;
t15 = -qJD(2) * pkin(3) + t77;
t12 = t36 * qJD(2) + t24;
t11 = t62 * t12;
t6 = t18 * qJD(5) - t28 * t63;
t5 = t74 * qJD(5) + t28 * t60;
t17 = [0, 0, -t61 * t123, -t64 * t123, -t21 * t28 + t22 * t29 + t24 * t32 + t129, t28 * qJD(2), t29 * qJD(2), t15 * t28 + t16 * t29 + t19 * t32 + t129, 0, 0, 0, 0, 0, -t6 * qJD(5) + (t63 * t101 + t29 * t60) * qJD(2), -t5 * qJD(5) + (-t60 * t101 + t29 * t63) * qJD(2), 0, 0, 0, 0, 0 (-t76 * qJD(6) + t29 * t62 - t5 * t59) * t49 + t75 * t86 + t6 * t37 - t74 * t26 -(t75 * qJD(6) + t29 * t59 + t5 * t62) * t49 - t76 * t86 + t6 * t39 - t74 * t25; 0, 0, 0, 0, t21 * t27 - t22 * t30 + (-t23 * t57 + t24 * t55) * pkin(2), 0 (0.2e1 * qJD(4) - t30) * qJD(2) + t24, -t15 * t27 + t97 * t16 + t19 * t50 + t23 * t88, -0.2e1 * t60 * t86, 0.2e1 * t112 * t96, -t116, -t115, 0, t78 * t104 + t68 * t60, -t78 * t105 + t68 * t63, t62 * t118 + t69 * t39 (t37 * t62 + t39 * t59) * t105 + (-t128 - t26 * t62 + (t37 * t59 - t39 * t62) * qJD(6)) * t63, t114 - t133, -t120 + (-t94 - t117) * qJD(5) + t72, t84 * t104 (-t34 * t103 - t113 * t62) * t49 + (-t87 * t102 + t37 * t106 + t67 * t59 + t11) * t60 + (t7 * t102 - t48 * t26 + t27 * t37 + t132 + (-t48 * t125 + (-t48 * t121 + t62 * t34) * qJD(2) + t1) * qJD(5)) * t63 (-t34 * t102 + t113 * t59) * t49 + (t39 * t106 + (t87 * qJD(6) - t12) * t59 + t67 * t62) * t60 + (-t7 * t103 - t48 * t25 + t27 * t39 + t131 + (-t48 * t119 - (t62 * t60 * t48 + t59 * t34) * qJD(2) - t2) * qJD(5)) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, t116, 0, 0, 0, 0, 0, t120 + (-t94 + t117) * qJD(5) + t72, t114 + t133; 0, 0, 0, 0, 0, 0, -t66, t135 * qJD(2), 0, 0, 0, 0, 0, t111 * t60, t111 * t63, 0, 0, 0, 0, 0, -t63 * t26 - t80 * t119 + (-t84 * t63 * t59 + t37 * t60) * qJD(5), -t118 + t80 * t125 + (-t63 * t119 + (t39 - t89) * t60) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, t63 * t66 * t60, -t112 * t66, 0, 0, 0, -t85 * t63, t85 * t60, t39 * t119 + t128 (t25 - t127) * t62 + (-t26 - t126) * t59, t92 + (t60 * t119 + (-t39 + t100) * t63) * qJD(2), -t49 * t103 + (-t49 * t121 + (t37 + t98) * t63) * qJD(2), -t49 * t107, -pkin(5) * t26 - t131 - (t62 * t42 - t59 * t9) * t49 - t10 * t37 + (-pkin(9) * t119 + t7 * t59) * qJD(6) + (-t1 * t63 + t70 * t59) * qJD(2), -pkin(5) * t25 + t132 + (t59 * t42 + t62 * t9) * t49 - t10 * t39 + (pkin(9) * t125 + t7 * t62) * qJD(6) + (t2 * t63 + t70 * t62) * qJD(2); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t37, -t37 ^ 2 + t39 ^ 2, t25 + t127, t126 - t26, t86, t134 * t2 - t59 * t3 - t7 * t39 + t11, t134 * t1 - t59 * t12 - t62 * t3 + t7 * t37;];
tauc_reg  = t17;
