% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% tauc_reg [5x18]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:20
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:19:59
% EndTime: 2022-01-20 10:20:02
% DurationCPUTime: 0.54s
% Computational Cost: add. (856->145), mult. (1680->189), div. (0->0), fcn. (927->6), ass. (0->105)
t65 = sin(qJ(4));
t67 = cos(qJ(4));
t60 = qJD(1) + qJD(2);
t100 = pkin(1) * qJD(1);
t68 = cos(qJ(2));
t87 = t68 * t100;
t44 = t60 * pkin(2) + t87;
t64 = cos(pkin(8));
t66 = sin(qJ(2));
t89 = t66 * t100;
t49 = t64 * t89;
t63 = sin(pkin(8));
t23 = t63 * t44 + t49;
t18 = t60 * pkin(7) + t23;
t98 = qJ(5) * t60;
t82 = t18 + t98;
t72 = t82 * t67;
t8 = t65 * qJD(3) + t72;
t56 = t67 * qJD(3);
t7 = -t82 * t65 + t56;
t97 = qJD(4) * pkin(4);
t6 = t7 + t97;
t117 = t6 - t7;
t61 = t65 ^ 2;
t116 = pkin(4) * t61;
t115 = t64 * pkin(2);
t114 = t67 * pkin(4);
t48 = t63 * t89;
t22 = t64 * t44 - t48;
t17 = -t60 * pkin(3) - t22;
t113 = t17 * t60;
t59 = t60 ^ 2;
t112 = t59 * t67;
t111 = t60 * t65;
t110 = t60 * t67;
t109 = t63 * t66;
t108 = t64 * t66;
t69 = qJD(4) ^ 2;
t107 = t69 * t65;
t85 = -pkin(3) - t114;
t11 = t85 * t60 + qJD(5) - t22;
t99 = pkin(1) * qJD(2);
t37 = (t63 * t68 + t108) * t99;
t31 = qJD(1) * t37;
t92 = t65 * qJD(4);
t86 = t60 * t92;
t16 = pkin(4) * t86 + t31;
t91 = t67 * qJD(4);
t106 = t11 * t91 + t16 * t65;
t105 = t17 * t91 + t31 * t65;
t36 = t63 * t87 + t49;
t38 = t64 * t87 - t48;
t104 = t36 * t110 + t38 * t92;
t54 = t68 * pkin(1) + pkin(2);
t103 = pkin(1) * t108 + t63 * t54;
t62 = t67 ^ 2;
t102 = -t61 - t62;
t101 = t61 - t62;
t35 = pkin(7) + t103;
t96 = -qJ(5) - t35;
t52 = t63 * pkin(2) + pkin(7);
t95 = -qJ(5) - t52;
t94 = qJD(5) * t60;
t90 = 0.2e1 * qJD(4) * t60;
t88 = pkin(4) * t92;
t14 = t18 * t92;
t39 = (t64 * t68 - t109) * t99;
t32 = qJD(1) * t39;
t78 = qJD(4) * qJD(3) + t32;
t2 = -qJ(5) * t86 - t14 + (t78 + t94) * t67;
t3 = (-t32 - t94) * t65 - t8 * qJD(4);
t84 = t2 * t67 - t3 * t65;
t83 = -pkin(1) * t109 + t64 * t54;
t81 = t67 * t90;
t80 = qJD(4) * t96;
t79 = qJD(4) * t95;
t34 = -pkin(3) - t83;
t77 = -t6 * t67 - t65 * t8;
t76 = t6 * t65 - t67 * t8;
t75 = (-qJD(2) + t60) * t100;
t74 = (-qJD(1) - t60) * t99;
t73 = t35 * t69 + t37 * t60;
t71 = qJD(4) * (t34 * t60 - t39);
t70 = (-qJD(5) - t11) * t60 - t78;
t58 = t67 * qJ(5);
t57 = t69 * t67;
t55 = t67 * qJD(5);
t53 = -pkin(3) - t115;
t45 = t85 - t115;
t43 = t65 * t81;
t42 = t67 * t52 + t58;
t41 = t95 * t65;
t33 = t101 * t90;
t30 = -t65 * qJD(5) + t67 * t79;
t29 = t65 * t79 + t55;
t28 = t34 - t114;
t27 = t38 * t91;
t25 = t37 + t88;
t20 = t67 * t35 + t58;
t19 = t96 * t65;
t12 = t17 * t92;
t9 = t11 * t92;
t5 = (-qJD(5) - t39) * t65 + t67 * t80;
t4 = t67 * t39 + t65 * t80 + t55;
t1 = [0, 0, 0, 0, t66 * t74, t68 * t74, t32 * t103 - t22 * t37 + t23 * t39 - t31 * t83, t43, -t33, t57, -t107, 0, t12 + t65 * t71 + (-t31 - t73) * t67, t73 * t65 + t67 * t71 + t105, t9 + (-t25 * t60 - t16) * t67 + (t28 * t111 + t5) * qJD(4), t25 * t111 + (t28 * t110 - t4) * qJD(4) + t106, (t4 * t67 - t5 * t65) * t60 + ((-t19 * t67 - t20 * t65) * t60 + t77) * qJD(4) + t84, t11 * t25 + t16 * t28 + t3 * t19 + t2 * t20 + t8 * t4 + t6 * t5; 0, 0, 0, 0, t66 * t75, t68 * t75, t22 * t36 - t23 * t38 + (-t31 * t64 + t32 * t63) * pkin(2), t43, -t33, t57, -t107, 0, t53 * t86 + t12 + (-t52 * t69 - t31) * t67 + t104, t52 * t107 + t27 + (-t36 * t65 + t53 * t91) * t60 + t105, -t16 * t67 + t9 + (t30 + (t45 - t114) * t111) * qJD(4) + t104, -t36 * t111 + t27 + (-t29 + (t45 * t67 + t116) * t60) * qJD(4) + t106, t77 * qJD(4) + (t29 * t67 - t30 * t65 + t102 * t38 + (-t41 * t67 - t42 * t65) * qJD(4)) * t60 + t84, t16 * t45 + t2 * t42 + t8 * t29 + t3 * t41 + t6 * t30 + t76 * t38 + (-t36 + t88) * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t57, -t107, -t57, 0, -t76 * qJD(4) + t2 * t65 + t3 * t67; 0, 0, 0, 0, 0, 0, 0, -t65 * t112, t101 * t59, 0, 0, 0, (-t32 - t113) * t65, t14 + (-t65 * t18 + t56) * qJD(4) + (-t78 - t113) * t67, (t8 - t72) * qJD(4) + (pkin(4) * t112 + t70) * t65, -t59 * t116 + t14 + (t65 * t98 + t7) * qJD(4) + t70 * t67, (-t97 + t117) * t110, t117 * t8 + (-t11 * t111 + t3) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t86, t81, t102 * t59, t76 * t60 + t16;];
tauc_reg = t1;
