% Calculate inertial parameters regressor of coriolis joint torque vector for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% tauc_reg [5x(5*10)]
%   inertial parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:13
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RPPRP1_coriolisvecJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_coriolisvecJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:12:45
% EndTime: 2022-01-23 09:12:47
% DurationCPUTime: 0.65s
% Computational Cost: add. (889->137), mult. (2126->211), div. (0->0), fcn. (1282->6), ass. (0->105)
t63 = sin(pkin(8));
t65 = cos(pkin(8));
t40 = -cos(pkin(7)) * pkin(1) - t65 * pkin(3) - t63 * pkin(6) - pkin(2);
t56 = sin(pkin(7)) * pkin(1) + qJ(3);
t68 = cos(qJ(4));
t44 = t68 * t65 * t56;
t67 = sin(qJ(4));
t19 = t67 * t40 + t44;
t30 = t40 * qJD(1) + qJD(3);
t50 = qJD(1) * t56;
t34 = t63 * qJD(2) + t65 * t50;
t11 = t67 * t30 + t68 * t34;
t72 = t11 * qJD(4);
t98 = qJD(1) * qJD(3);
t84 = t65 * t98;
t5 = -t67 * t84 - t72;
t10 = t68 * t30 - t67 * t34;
t75 = t10 * t67 - t11 * t68;
t102 = qJD(4) * t68;
t103 = qJD(4) * t67;
t94 = -t30 * t102 + t34 * t103 - t68 * t84;
t119 = -t75 * qJD(4) + t5 * t68 - t67 * t94;
t101 = qJD(5) * t67;
t107 = qJD(1) * t63;
t85 = qJ(5) * t102;
t1 = (-t85 - t101) * t107 - t94;
t89 = t63 * t103;
t81 = qJD(1) * t89;
t48 = qJ(5) * t81;
t113 = t63 * t68;
t105 = qJD(3) * t67;
t87 = t65 * t105;
t71 = -qJD(5) * t113 - t87;
t2 = t71 * qJD(1) + t48 - t72;
t100 = t65 * qJD(1);
t55 = -qJD(4) + t100;
t86 = qJ(5) * t107;
t8 = -t68 * t86 + t10;
t3 = -t55 * pkin(4) + t8;
t9 = -t67 * t86 + t11;
t77 = t3 * t67 - t68 * t9;
t118 = -t77 * qJD(4) + t1 * t67 + t2 * t68;
t117 = t3 - t8;
t58 = t65 * qJD(2);
t33 = t63 * t50 - t58;
t116 = t33 * t63;
t115 = t56 * t67;
t59 = t63 ^ 2;
t69 = qJD(1) ^ 2;
t114 = t59 * t69;
t104 = qJD(3) * t68;
t112 = t40 * t102 + t65 * t104;
t97 = qJD(1) * qJD(4);
t83 = t68 * t97;
t80 = t63 * t83;
t38 = pkin(4) * t80 + t63 * t98;
t60 = t65 ^ 2;
t111 = t59 + t60;
t61 = t67 ^ 2;
t62 = t68 ^ 2;
t110 = t61 - t62;
t109 = qJ(5) * t63;
t108 = qJD(1) * t59;
t106 = qJD(1) * t67;
t20 = qJD(5) - t58 + (pkin(4) * t106 + t50) * t63;
t99 = qJD(5) + t20;
t96 = t67 * t114;
t95 = t65 * t115;
t93 = t68 * t109;
t92 = t63 * t106;
t91 = t68 * t107;
t90 = t56 * t103;
t88 = t63 * t102;
t82 = t55 * t89;
t78 = t3 * t68 + t67 * t9;
t74 = -t34 * t65 - t116;
t73 = t59 * t67 * t83;
t70 = -t55 ^ 2 - t114;
t49 = t68 * t96;
t47 = t65 * t81;
t46 = -0.2e1 * t73;
t45 = 0.2e1 * t73;
t43 = t59 * t56 * t98;
t42 = (pkin(4) * t102 + qJD(3)) * t63;
t41 = t55 * t91;
t39 = t110 * t114;
t37 = (pkin(4) * t67 + t56) * t63;
t36 = t68 * t40;
t31 = 0.2e1 * t110 * t59 * t97;
t28 = -t41 - t80;
t27 = (-qJD(4) - t55) * t92;
t24 = (t55 - t100) * t88;
t23 = (t55 + t100) * t88;
t22 = t47 + t82;
t21 = t47 - t82;
t18 = t36 - t95;
t17 = t70 * t68;
t16 = t70 * t67;
t15 = -t67 * t109 + t19;
t14 = -t19 * qJD(4) - t87;
t13 = -t65 * t90 + t112;
t12 = -t93 + t36 + (-pkin(4) - t115) * t65;
t7 = (-t44 + (-t40 + t109) * t67) * qJD(4) + t71;
t6 = -t63 * t101 + (-t93 - t95) * qJD(4) + t112;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t111 * t98, t43 + (t60 * t50 - t74) * qJD(3), t46, t31, t22, t45, t23, 0, t33 * t88 - t14 * t55 - t5 * t65 + (t56 * t102 + 0.2e1 * t105) * t108, -t33 * t89 + t13 * t55 - t94 * t65 + (-t90 + 0.2e1 * t104) * t108, ((-t13 * t67 - t14 * t68 + (t18 * t67 - t19 * t68) * qJD(4)) * qJD(1) - t119) * t63, qJD(3) * t116 + t10 * t14 + t11 * t13 + t5 * t18 - t19 * t94 + t43, t46, t31, t22, t45, t23, 0, -t2 * t65 - t7 * t55 + (t20 * t102 + t38 * t67 + (t37 * t102 + t42 * t67) * qJD(1)) * t63, t1 * t65 + t6 * t55 + (-t20 * t103 + t38 * t68 + (-t37 * t103 + t42 * t68) * qJD(1)) * t63, ((-t6 * t67 - t68 * t7 + (t12 * t67 - t15 * t68) * qJD(4)) * qJD(1) - t118) * t63, t1 * t15 + t2 * t12 + t20 * t42 + t3 * t7 + t38 * t37 + t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t21, 0, (-t84 - t94 * t68 - t5 * t67 + (-t10 * t68 - t11 * t67) * qJD(4)) * t63, 0, 0, 0, 0, 0, 0, t24, t21, 0, -t38 * t65 + (-t78 * qJD(4) + t1 * t68 - t2 * t67) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111 * t69, t74 * qJD(1), 0, 0, 0, 0, 0, 0, t16, t17, 0, (t75 * t65 - t116) * qJD(1) + t119, 0, 0, 0, 0, 0, 0, t16, t17, 0, (-t20 * t63 + t77 * t65) * qJD(1) + t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t39, t27, -t49, t28, 0, -t11 * t55 - t72 + (-t33 * t113 - t87) * qJD(1), -t10 * t55 + t33 * t92 + t94, 0, 0, t49, -t39, t27, -t49, t28, 0, -t9 * t55 + t48 + (-qJD(4) * t30 - t84) * t67 + (-pkin(4) * t96 - qJD(4) * t34 - t99 * t107) * t68, -t62 * pkin(4) * t114 - t8 * t55 + (t99 * t67 + t85) * t107 + t94, (pkin(4) * qJD(4) - t117) * t92, t117 * t9 + (-t20 * t91 + t2) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41 + t80, (-qJD(4) + t55) * t92, (-t61 - t62) * t114, t78 * t107 + t38;];
tauc_reg = t4;
