% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:55
% DurationCPUTime: 0.77s
% Computational Cost: add. (1512->143), mult. (2975->215), div. (0->0), fcn. (1780->8), ass. (0->113)
t70 = sin(pkin(8));
t109 = pkin(1) * qJD(1);
t74 = sin(qJ(2));
t99 = t74 * t109;
t56 = t70 * t99;
t72 = cos(pkin(8));
t76 = cos(qJ(2));
t98 = t76 * t109;
t44 = t72 * t98 - t56;
t101 = qJD(4) - t44;
t71 = cos(pkin(9));
t75 = cos(qJ(5));
t117 = t71 * t75;
t39 = t44 * qJD(2);
t66 = qJD(1) + qJD(2);
t32 = t66 * qJD(4) + t39;
t116 = t72 * t74;
t82 = t70 * t76 + t116;
t80 = pkin(1) * t82;
t43 = qJD(2) * t80;
t38 = qJD(1) * t43;
t69 = sin(pkin(9));
t81 = -t71 * pkin(4) - t69 * pkin(7) - pkin(3);
t51 = t66 * pkin(2) + t98;
t34 = t72 * t51 - t56;
t87 = qJD(4) - t34;
t12 = t81 * t66 + t87;
t35 = t70 * t51 + t72 * t99;
t30 = t66 * qJ(4) + t35;
t19 = t69 * qJD(3) + t71 * t30;
t73 = sin(qJ(5));
t7 = t75 * t12 - t73 * t19;
t2 = t7 * qJD(5) + t32 * t117 + t73 * t38;
t131 = t2 * t73;
t130 = t72 * pkin(2);
t64 = t69 ^ 2;
t121 = t64 * t75;
t129 = t32 * t121 + t2 * t71;
t18 = -t71 * qJD(3) + t69 * t30;
t128 = t18 * t69;
t127 = t32 * t71;
t108 = pkin(1) * qJD(2);
t58 = t70 * t74 * pkin(1);
t45 = t72 * t76 * t108 - qJD(2) * t58;
t40 = qJD(4) + t45;
t126 = t40 * t66;
t42 = qJD(1) * t80;
t125 = t42 * t66;
t124 = t43 * t66;
t28 = t64 * t32;
t63 = t66 ^ 2;
t123 = t64 * t63;
t122 = t64 * t66;
t65 = t71 ^ 2;
t29 = t65 * t32;
t120 = t66 * t69;
t119 = t71 * t66;
t118 = t71 * t73;
t102 = t71 * qJD(4);
t46 = t81 - t130;
t60 = t70 * pkin(2) + qJ(4);
t24 = t60 * t117 + t73 * t46;
t105 = qJD(5) * t24;
t115 = t73 * t102 - t44 * t118 + t75 * t42 + t105;
t23 = -t60 * t118 + t75 * t46;
t106 = qJD(5) * t23;
t114 = t75 * t102 - t44 * t117 - t73 * t42 + t106;
t113 = t28 + t29;
t61 = t76 * pkin(1) + pkin(2);
t112 = pkin(1) * t116 + t70 * t61;
t111 = t64 + t65;
t110 = t73 ^ 2 - t75 ^ 2;
t8 = t73 * t12 + t75 * t19;
t107 = qJD(5) * t8;
t104 = qJD(5) * t73;
t103 = qJD(5) * t75;
t100 = t18 * t120;
t97 = qJD(5) * t122;
t96 = t69 * t104;
t95 = t69 * t103;
t94 = t72 * t61 - t58;
t93 = t63 * t73 * t121;
t55 = -qJD(5) + t119;
t91 = t55 * t96;
t3 = -t32 * t118 + t75 * t38 - t107;
t90 = t18 * t95 + t73 * t28 - t3 * t71;
t89 = (-qJD(5) - t55) * t120;
t88 = t7 * t73 - t75 * t8;
t86 = (-qJD(2) + t66) * t109;
t85 = (-qJD(1) - t66) * t108;
t84 = t75 * t73 * t97;
t83 = t19 * t71 + t128;
t31 = t81 - t94;
t41 = qJ(4) + t112;
t10 = t41 * t117 + t73 * t31;
t9 = -t41 * t118 + t75 * t31;
t79 = -t131 + (-t3 - t107) * t75;
t78 = -t55 ^ 2 - t123;
t77 = t82 * qJD(1) * t108;
t50 = t96 * t119;
t49 = -0.2e1 * t84;
t48 = 0.2e1 * t84;
t36 = t38 * t69;
t33 = 0.2e1 * t110 * t97;
t27 = -t66 * pkin(3) + t87;
t26 = (t55 + t119) * t95;
t25 = t50 + t91;
t20 = t60 * t28;
t11 = t41 * t28;
t6 = t7 * t96;
t5 = -t10 * qJD(5) - t40 * t118 + t75 * t43;
t4 = t9 * qJD(5) + t40 * t117 + t73 * t43;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t85, t76 * t85, 0, 0, 0, 0, 0, 0, 0, 0, -t77 - t124, -t45 * t66 - t39, 0, t39 * t112 - t34 * t43 + t35 * t45 - t38 * t94, 0, 0, 0, 0, 0, 0, (-t38 - t124) * t71, t43 * t120 + t36, t111 * t126 + t113, t41 * t29 + t11 + t38 * (-pkin(3) - t94) + t27 * t43 + t83 * t40, t49, t33, t25, t48, t26, 0, -t5 * t55 + (t41 * t103 + t40 * t73) * t122 + t90, t121 * t126 + t4 * t55 + (-t41 * t122 - t128) * t104 + t129, t6 + ((-t4 * t73 - t5 * t75 + (-t10 * t75 + t73 * t9) * qJD(5)) * t66 + t79) * t69, t2 * t10 + t40 * t128 + t3 * t9 + t8 * t4 + t7 * t5 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74 * t86, t76 * t86, 0, 0, 0, 0, 0, 0, 0, 0, -t77 + t125, t44 * t66 - t39, 0, t34 * t42 - t35 * t44 + (-t38 * t72 + t39 * t70) * pkin(2), 0, 0, 0, 0, 0, 0, (-t38 + t125) * t71, -t42 * t120 + t36, t101 * t66 * t111 + t113, t60 * t29 + t20 + t38 * (-pkin(3) - t130) - t27 * t42 + t101 * t83, t49, t33, t25, t48, t26, 0, t115 * t55 + (t101 * t73 + t60 * t103) * t122 + t90, -t18 * t96 + t114 * t55 + (t101 * t75 - t60 * t104) * t122 + t129, t6 + (((-t105 + t115) * t75 + (t106 - t114) * t73) * t66 + t79) * t69, t101 * t128 + t114 * t8 - t115 * t7 + t2 * t24 + t3 * t23 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t55 - t119) * t95, t50 - t91, 0, (t2 * t75 - t3 * t73 - t127 + (-t7 * t75 - t73 * t8) * qJD(5)) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t63, -t83 * t66 + t38, 0, 0, 0, 0, 0, 0, t78 * t73, t78 * t75, 0, t131 + t3 * t75 - t88 * qJD(5) + (t88 * t71 - t128) * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t93, -t110 * t123, t73 * t89, -t93, t75 * t89, 0, -t75 * t100 - t8 * t55 + t3, -t7 * t55 + (-qJD(5) * t12 - t127) * t75 + (qJD(5) * t19 + t100 - t38) * t73, 0, 0;];
tauc_reg = t1;
