% Calculate minimal parameter regressor of coriolis joint torque vector for
% S4RRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% tauc_reg [4x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S4RRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:23:18
% EndTime: 2019-12-31 17:23:21
% DurationCPUTime: 0.62s
% Computational Cost: add. (603->107), mult. (1101->175), div. (0->0), fcn. (698->6), ass. (0->97)
t118 = pkin(6) + pkin(7);
t62 = sin(qJ(4));
t63 = sin(qJ(3));
t65 = cos(qJ(4));
t66 = cos(qJ(3));
t36 = t62 * t66 + t65 * t63;
t59 = qJD(1) + qJD(2);
t28 = t36 * t59;
t58 = qJD(3) + qJD(4);
t117 = qJD(4) - t58;
t67 = cos(qJ(2));
t115 = t67 * pkin(1);
t64 = sin(qJ(2));
t52 = t64 * pkin(1) + pkin(6);
t114 = -pkin(7) - t52;
t16 = t58 * t36;
t54 = -t66 * pkin(3) - pkin(2);
t99 = pkin(1) * qJD(1);
t87 = t67 * t99;
t29 = t54 * t59 - t87;
t98 = pkin(1) * qJD(2);
t83 = qJD(1) * t98;
t78 = t64 * t83;
t97 = qJD(3) * t63;
t86 = t59 * t97;
t30 = pkin(3) * t86 + t78;
t103 = t65 * t66;
t108 = t62 * t63;
t35 = -t103 + t108;
t113 = t29 * t16 + t30 * t35;
t15 = t58 * t35;
t112 = -t29 * t15 + t30 * t36;
t90 = t59 * t103;
t91 = t59 * t108;
t26 = -t90 + t91;
t111 = t28 * t26;
t110 = t58 * t67;
t109 = t59 * t63;
t106 = t64 * t66;
t81 = t118 * t59 + t64 * t99;
t24 = t81 * t66;
t105 = t65 * t24;
t68 = qJD(3) ^ 2;
t102 = t68 * t63;
t55 = t68 * t66;
t42 = -t59 * pkin(2) - t87;
t96 = qJD(3) * t66;
t101 = t42 * t96 + t63 * t78;
t100 = t63 ^ 2 - t66 ^ 2;
t95 = qJD(3) * t67;
t94 = -qJD(1) - t59;
t93 = -qJD(2) + t59;
t92 = pkin(3) * t109;
t89 = t67 * t98;
t88 = pkin(3) * t97;
t85 = t59 * t96;
t84 = t63 * t95;
t23 = t81 * t63;
t20 = qJD(3) * pkin(3) - t23;
t82 = -pkin(3) * t58 - t20;
t80 = qJD(3) * t118;
t79 = qJD(3) * t114;
t77 = t67 * t83;
t76 = t93 * t99;
t75 = t94 * t98;
t73 = qJD(3) * t81;
t11 = -t63 * t73 + t66 * t77;
t12 = -t63 * t77 - t66 * t73;
t72 = -t62 * t11 + t65 * t12 - t29 * t28;
t71 = -t64 * t109 + t66 * t95;
t70 = -t42 * t59 - t77;
t8 = qJD(4) * t90 - t58 * t91 + t65 * t85;
t69 = t29 * t26 + (t117 * t24 - t12) * t62;
t57 = t59 ^ 2;
t56 = t66 * pkin(7);
t53 = -pkin(2) - t115;
t51 = t66 * pkin(6) + t56;
t50 = t118 * t63;
t45 = t54 - t115;
t40 = 0.2e1 * t63 * t85;
t39 = t64 * t98 + t88;
t38 = t66 * t80;
t37 = t63 * t80;
t34 = t66 * t52 + t56;
t33 = t114 * t63;
t31 = t42 * t97;
t25 = -0.2e1 * t100 * t59 * qJD(3);
t22 = -t63 * t89 + t66 * t79;
t21 = t63 * t79 + t66 * t89;
t14 = t16 * t58;
t13 = t15 * t58;
t9 = t16 * t59;
t5 = -t26 ^ 2 + t28 ^ 2;
t3 = t26 * t58 + t8;
t2 = -t28 * t15 + t8 * t36;
t1 = t15 * t26 - t28 * t16 - t8 * t35 - t36 * t9;
t4 = [0, 0, 0, 0, t64 * t75, t67 * t75, t40, t25, t55, -t102, 0, t53 * t86 - t52 * t55 + t31 + (t94 * t106 - t84) * t98, t52 * t102 + t53 * t85 - t71 * t98 + t101, t2, t1, -t13, -t14, 0, t39 * t26 + t45 * t9 + (-t62 * t21 + t65 * t22 + (-t33 * t62 - t34 * t65) * qJD(4)) * t58 + t113, t39 * t28 + t45 * t8 - (t65 * t21 + t62 * t22 + (t33 * t65 - t34 * t62) * qJD(4)) * t58 + t112; 0, 0, 0, 0, t64 * t76, t67 * t76, t40, t25, t55, -t102, 0, -pkin(2) * t86 - pkin(6) * t55 + t31 + (t93 * t106 + t84) * t99, -pkin(2) * t85 + pkin(6) * t102 + t71 * t99 + t101, t2, t1, -t13, -t14, 0, t26 * t88 + t54 * t9 + (t62 * t37 - t65 * t38 + (t50 * t62 - t51 * t65) * qJD(4)) * t58 + (t36 * t110 - t64 * t26) * t99 + t113, t28 * t88 + t54 * t8 - (-t65 * t37 - t62 * t38 + (-t50 * t65 - t51 * t62) * qJD(4)) * t58 + (-t35 * t110 - t64 * t28) * t99 + t112; 0, 0, 0, 0, 0, 0, -t63 * t57 * t66, t100 * t57, 0, 0, 0, t70 * t63, t70 * t66, t111, t5, t3, 0, 0, -t26 * t92 - (t62 * t23 - t105) * t58 + (t82 * t62 - t105) * qJD(4) + t72, -t28 * t92 + (t82 * qJD(4) - t23 * t58 - t11) * t65 + t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111, t5, t3, 0, 0, t72 + t117 * (-t62 * t20 - t105), (-t117 * t20 - t11) * t65 + t69;];
tauc_reg = t4;
