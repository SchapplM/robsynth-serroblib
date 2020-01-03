% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRPRP2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:43
% EndTime: 2019-12-31 19:49:44
% DurationCPUTime: 0.43s
% Computational Cost: add. (818->122), mult. (1642->157), div. (0->0), fcn. (889->6), ass. (0->92)
t58 = sin(pkin(8));
t61 = sin(qJ(2));
t103 = t58 * t61;
t59 = cos(pkin(8));
t63 = cos(qJ(2));
t91 = pkin(1) * qJD(2);
t35 = (t59 * t63 - t103) * t91;
t28 = qJD(1) * t35;
t112 = qJD(3) * qJD(4) + t28;
t55 = qJD(1) + qJD(2);
t92 = pkin(1) * qJD(1);
t82 = t63 * t92;
t41 = t55 * pkin(2) + t82;
t83 = t61 * t92;
t44 = t59 * t83;
t22 = t58 * t41 + t44;
t18 = t55 * pkin(7) + t22;
t60 = sin(qJ(4));
t101 = t60 * t18;
t62 = cos(qJ(4));
t11 = t62 * qJD(3) - t101;
t111 = qJD(5) - t11;
t97 = t112 * t62;
t2 = (qJD(5) - t101) * qJD(4) + t97;
t89 = qJD(4) * t62;
t4 = t112 * t60 + t18 * t89;
t110 = t2 * t62 + t4 * t60;
t56 = t60 ^ 2;
t57 = t62 ^ 2;
t109 = (t56 + t57) * t55;
t9 = -qJD(4) * pkin(4) + t111;
t12 = t60 * qJD(3) + t62 * t18;
t86 = qJD(4) * qJ(5);
t10 = t12 + t86;
t108 = t59 * pkin(2);
t102 = t59 * t61;
t51 = t63 * pkin(1) + pkin(2);
t95 = pkin(1) * t102 + t58 * t51;
t31 = pkin(7) + t95;
t64 = qJD(4) ^ 2;
t107 = t31 * t64;
t47 = t58 * pkin(2) + pkin(7);
t106 = t47 * t64;
t105 = t55 * t60;
t104 = t55 * t62;
t100 = t64 * t60;
t43 = t58 * t83;
t21 = t59 * t41 - t43;
t17 = -t55 * pkin(3) - t21;
t33 = (t58 * t63 + t102) * t91;
t27 = qJD(1) * t33;
t99 = t17 * t89 + t27 * t60;
t32 = t58 * t82 + t44;
t34 = t59 * t82 - t43;
t90 = qJD(4) * t60;
t98 = t32 * t104 + t34 * t90;
t87 = t60 * qJD(5);
t36 = pkin(4) * t90 - t62 * t86 - t87;
t96 = t32 - t36;
t94 = t56 - t57;
t88 = t10 * qJD(4);
t54 = t55 ^ 2;
t84 = t60 * t54 * t62;
t81 = t55 * t90;
t72 = pkin(4) * t60 - qJ(5) * t62;
t6 = t27 + (t72 * qJD(4) - t87) * t55;
t80 = -t6 - t106;
t79 = t11 + t101;
t68 = -t62 * pkin(4) - t60 * qJ(5) - pkin(3);
t77 = -pkin(1) * t103 + t59 * t51;
t19 = t68 - t77;
t78 = t19 * t55 - t35;
t75 = (-qJD(2) + t55) * t92;
t74 = (-qJD(1) - t55) * t91;
t73 = t10 * t62 + t60 * t9;
t71 = t33 * t55 + t107;
t70 = qJD(4) * ((-pkin(3) - t77) * t55 - t35);
t69 = t12 * qJD(4) - t4;
t13 = t33 + t36;
t67 = -t13 * t55 - t107 - t6;
t66 = -t60 * t88 + t9 * t89 + t110;
t65 = (-t10 * t60 + t62 * t9) * qJD(4) + t110;
t53 = t64 * t62;
t48 = -pkin(3) - t108;
t40 = 0.2e1 * t62 * t81;
t39 = t68 - t108;
t37 = t72 * t55;
t29 = -0.2e1 * t94 * t55 * qJD(4);
t14 = t17 * t90;
t8 = t68 * t55 - t21;
t5 = t8 * t90;
t1 = [0, 0, 0, 0, t61 * t74, t63 * t74, -t21 * t33 + t22 * t35 - t27 * t77 + t28 * t95, t40, t29, t53, -t100, 0, t14 + t60 * t70 + (-t27 - t71) * t62, t71 * t60 + t62 * t70 + t99, t67 * t62 + t78 * t90 + t5, t35 * t109 + t66, t67 * t60 + (-t78 - t8) * t89, t8 * t13 + t6 * t19 + t65 * t31 + t73 * t35; 0, 0, 0, 0, t61 * t75, t63 * t75, t21 * t32 - t22 * t34 + (-t27 * t59 + t28 * t58) * pkin(2), t40, t29, t53, -t100, 0, t48 * t81 + t14 + (-t27 - t106) * t62 + t98, (-t32 * t55 + t106) * t60 + (t48 * t55 + t34) * t89 + t99, t39 * t81 + t5 + (-t36 * t55 + t80) * t62 + t98, -t34 * t109 + t66, (-t39 * t55 - t34 - t8) * t89 + (t96 * t55 + t80) * t60, -t73 * t34 + t6 * t39 + t65 * t47 - t96 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, -t53, -t100, 0, t53, t73 * qJD(4) + t2 * t60 - t4 * t62; 0, 0, 0, 0, 0, 0, 0, -t84, t94 * t54, 0, 0, 0, -t17 * t105 + t69, t79 * qJD(4) - t17 * t104 - t97, (t37 * t62 - t60 * t8) * t55 + t69, 0, (t37 * t60 + t62 * t8) * t55 + (0.2e1 * qJD(5) - t79) * qJD(4) + t97, -t4 * pkin(4) + t2 * qJ(5) + t111 * t10 - t9 * t12 - t8 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, -t56 * t54 - t64, t8 * t105 + t4 - t88;];
tauc_reg = t1;
