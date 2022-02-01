% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:31
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRPR2_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:30:53
% EndTime: 2022-01-20 11:30:55
% DurationCPUTime: 0.73s
% Computational Cost: add. (2042->135), mult. (4308->195), div. (0->0), fcn. (2463->8), ass. (0->106)
t107 = pkin(2) * qJD(3);
t66 = sin(pkin(9));
t69 = sin(qJ(3));
t120 = t66 * t69;
t108 = pkin(1) * qJD(1);
t70 = sin(qJ(2));
t72 = cos(qJ(3));
t117 = t70 * t72;
t73 = cos(qJ(2));
t88 = -t69 * t73 - t117;
t44 = t88 * t108;
t118 = t69 * t70;
t87 = t72 * t73 - t118;
t45 = t87 * t108;
t67 = cos(pkin(9));
t113 = -t66 * t44 - t67 * t45 + (t67 * t72 - t120) * t107;
t63 = qJD(1) + qJD(2);
t61 = qJD(3) + t63;
t132 = t113 * t61;
t101 = t73 * t108;
t50 = pkin(2) * t63 + t101;
t102 = t70 * t108;
t98 = t69 * t102;
t33 = t72 * t50 - t98;
t30 = pkin(3) * t61 + t33;
t34 = t102 * t72 + t50 * t69;
t31 = t67 * t34;
t16 = t30 * t66 + t31;
t14 = pkin(8) * t61 + t16;
t68 = sin(qJ(5));
t71 = cos(qJ(5));
t10 = qJD(4) * t68 + t14 * t71;
t112 = (qJD(2) + qJD(3)) * t98;
t96 = qJD(2) * t101;
t20 = (qJD(3) * t50 + t96) * t72 - t112;
t106 = qJD(3) * t69;
t105 = qJD(3) * t72;
t76 = (qJD(2) * t88 - t105 * t70) * pkin(1);
t75 = qJD(1) * t76;
t21 = -t106 * t50 + t75;
t6 = t67 * t20 + t66 * t21;
t9 = qJD(4) * t71 - t14 * t68;
t2 = t9 * qJD(5) + t71 * t6;
t3 = -t10 * qJD(5) - t68 * t6;
t131 = t2 * t71 - t3 * t68 + (-t10 * t68 - t71 * t9) * qJD(5);
t5 = t66 * t20 - t21 * t67;
t119 = t67 * t69;
t114 = -t67 * t44 + t66 * t45 - (t66 * t72 + t119) * t107;
t99 = t114 * t61;
t130 = t99 - t5;
t92 = t10 * t71 - t68 * t9;
t59 = pkin(1) * t73 + pkin(2);
t27 = t59 * t105 + (qJD(2) * t87 - t106 * t70) * pkin(1);
t28 = -t106 * t59 + t76;
t7 = t27 * t66 - t28 * t67;
t128 = t7 * t61;
t8 = t27 * t67 + t28 * t66;
t127 = t8 * t61;
t104 = qJD(5) * t71;
t121 = t66 * t34;
t15 = t30 * t67 - t121;
t13 = -pkin(4) * t61 - t15;
t126 = t104 * t13 + t5 * t68;
t17 = t33 * t66 + t31;
t124 = t17 * t61;
t18 = t33 * t67 - t121;
t123 = t18 * t61;
t58 = pkin(2) * t72 + pkin(3);
t111 = pkin(2) * t119 + t66 * t58;
t40 = pkin(8) + t111;
t74 = qJD(5) ^ 2;
t122 = t40 * t74;
t116 = t74 * t68;
t95 = -pkin(1) * t118 + t72 * t59;
t43 = pkin(3) + t95;
t46 = pkin(1) * t117 + t59 * t69;
t115 = t43 * t66 + t67 * t46;
t64 = t68 ^ 2;
t65 = t71 ^ 2;
t110 = t64 - t65;
t109 = t64 + t65;
t60 = t61 ^ 2;
t103 = t68 * t60 * t71;
t100 = -t13 * t61 - t6;
t97 = t68 * t61 * t104;
t94 = (-qJD(2) + t63) * t108;
t93 = pkin(1) * qJD(2) * (-qJD(1) - t63);
t23 = pkin(8) + t115;
t91 = t23 * t74 + t128;
t56 = pkin(3) * t66 + pkin(8);
t90 = t56 * t74 - t124;
t89 = t43 * t67 - t46 * t66;
t22 = -pkin(4) - t89;
t86 = qJD(5) * (t22 * t61 - t8);
t85 = (-pkin(2) * t61 - t50) * qJD(3);
t57 = -pkin(3) * t67 - pkin(4);
t84 = qJD(5) * (t57 * t61 + t18);
t83 = -pkin(2) * t120 + t58 * t67;
t39 = -pkin(4) - t83;
t82 = qJD(5) * (t39 * t61 - t113);
t62 = t74 * t71;
t49 = -0.2e1 * t97;
t48 = 0.2e1 * t97;
t35 = -0.2e1 * t110 * t61 * qJD(5);
t11 = t13 * qJD(5) * t68;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t93, t73 * t93, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t61 + t21, -t27 * t61 - t20, 0, t20 * t46 + t21 * t95 + t27 * t34 + t28 * t33, 0, 0, 0, 0, 0, 0, -t5 - t128, -t6 - t127, 0, t115 * t6 - t15 * t7 + t16 * t8 - t5 * t89, t48, t35, t62, t49, -t116, 0, t11 + t68 * t86 + (-t5 - t91) * t71, t68 * t91 + t71 * t86 + t126, t109 * t127 + t131, t13 * t7 + t131 * t23 + t5 * t22 + t8 * t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70 * t94, t73 * t94, 0, 0, 0, 0, 0, 0, 0, 0, -t44 * t61 + t69 * t85 + t75, t45 * t61 + (t85 - t96) * t72 + t112, 0, -t33 * t44 - t34 * t45 + (t20 * t69 + t21 * t72 + (-t33 * t69 + t34 * t72) * qJD(3)) * pkin(2), 0, 0, 0, 0, 0, 0, t130, -t6 - t132, 0, t111 * t6 + t113 * t16 + t114 * t15 - t5 * t83, t48, t35, t62, t49, -t116, 0, t11 + t68 * t82 + (-t122 + t130) * t71, (t122 - t99) * t68 + t71 * t82 + t126, t109 * t132 + t131, t113 * t92 - t114 * t13 + t131 * t40 + t5 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t61 + t21, t33 * t61 - t20, 0, 0, 0, 0, 0, 0, 0, 0, -t5 + t124, -t6 + t123, 0, t15 * t17 - t16 * t18 + (-t5 * t67 + t6 * t66) * pkin(3), t48, t35, t62, t49, -t116, 0, t11 + t68 * t84 + (-t5 - t90) * t71, t68 * t90 + t71 * t84 + t126, -t109 * t123 + t131, -t13 * t17 + t131 * t56 - t18 * t92 + t5 * t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116, -t62, 0, qJD(5) * t92 + t2 * t68 + t3 * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t103, t110 * t60, 0, t103, 0, 0, t100 * t68, t100 * t71, 0, 0;];
tauc_reg = t1;
