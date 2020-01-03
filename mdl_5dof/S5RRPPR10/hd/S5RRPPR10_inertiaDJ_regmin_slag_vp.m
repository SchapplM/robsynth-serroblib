% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR10_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:44:40
% EndTime: 2019-12-31 19:44:43
% DurationCPUTime: 0.81s
% Computational Cost: add. (543->141), mult. (1460->270), div. (0->0), fcn. (1171->6), ass. (0->87)
t53 = sin(pkin(8));
t54 = cos(pkin(8));
t55 = sin(qJ(5));
t57 = cos(qJ(5));
t33 = t57 * t53 - t55 * t54;
t56 = sin(qJ(2));
t24 = t33 * t56;
t51 = t53 ^ 2;
t99 = (t54 ^ 2 + t51) * qJD(3);
t58 = cos(qJ(2));
t77 = t56 * qJD(2);
t98 = qJ(4) * t77 - t58 * qJD(4);
t69 = t53 * qJ(4) + pkin(2);
t34 = -t54 * pkin(3) - t69;
t85 = t56 * qJ(3);
t97 = qJD(2) * (-t34 * t58 + t85);
t96 = pkin(3) + pkin(4);
t95 = t53 * t56;
t94 = t53 * t58;
t26 = -t56 * qJD(3) + (pkin(2) * t56 - qJ(3) * t58) * qJD(2);
t93 = t54 * t26;
t92 = t54 * t56;
t91 = t54 * t58;
t88 = -pkin(7) + qJ(3);
t67 = -t58 * pkin(2) - t85;
t35 = -pkin(1) + t67;
t45 = pkin(6) * t91;
t21 = t53 * t35 + t45;
t87 = qJ(3) * t99;
t86 = qJ(4) * t54;
t84 = qJD(3) * t58;
t83 = qJD(5) * t55;
t82 = qJD(5) * t57;
t81 = qJD(5) * t58;
t80 = t51 * qJD(4);
t79 = t53 * qJD(3);
t78 = t53 * qJD(4);
t76 = t58 * qJD(2);
t44 = pkin(6) * t94;
t73 = -0.2e1 * pkin(1) * qJD(2);
t72 = pkin(6) * t77;
t49 = pkin(6) * t76;
t43 = t53 * t76;
t41 = t54 * t76;
t71 = t56 * t76;
t70 = -pkin(6) * t53 - pkin(3);
t20 = t54 * t35 - t44;
t18 = -t58 * qJ(4) + t21;
t68 = -qJD(4) * t92 + t49;
t66 = pkin(3) * t53 - t86;
t50 = t58 * pkin(3);
t11 = t58 * pkin(4) + t44 + t50 + (-pkin(7) * t56 - t35) * t54;
t13 = pkin(7) * t95 + t18;
t65 = t57 * t11 - t55 * t13;
t64 = t55 * t11 + t57 * t13;
t15 = t53 * t72 + t93;
t22 = t53 * t26;
t16 = -t54 * t72 + t22;
t63 = -t15 * t53 + t16 * t54;
t37 = t88 * t53;
t38 = t88 * t54;
t62 = t57 * t37 - t55 * t38;
t61 = t55 * t37 + t57 * t38;
t32 = t55 * t53 + t57 * t54;
t59 = -t96 * t53 + t86;
t27 = t32 * qJD(5);
t40 = t58 * t79;
t31 = 0.2e1 * t99;
t29 = t96 * t54 + t69;
t28 = t53 * t82 - t54 * t83;
t25 = t32 * t56;
t23 = (pkin(6) + t66) * t56;
t19 = -t20 + t50;
t17 = (-pkin(6) + t59) * t56;
t14 = t66 * t76 + t68;
t12 = t70 * t77 - t93;
t10 = t59 * t76 - t68;
t9 = t16 + t98;
t8 = qJD(5) * t24 + t32 * t76;
t7 = t27 * t56 + t55 * t41 - t57 * t43;
t6 = qJD(3) * t33 - qJD(5) * t61;
t5 = -qJD(3) * t32 - qJD(5) * t62;
t4 = t22 + (-pkin(6) * t92 + pkin(7) * t94) * qJD(2) + t98;
t3 = -t93 + (-pkin(7) * t91 + (-pkin(4) + t70) * t56) * qJD(2);
t2 = -qJD(5) * t64 + t57 * t3 - t55 * t4;
t1 = -qJD(5) * t65 - t55 * t3 - t57 * t4;
t30 = [0, 0, 0, 0.2e1 * t71, 0.2e1 * (-t56 ^ 2 + t58 ^ 2) * qJD(2), 0, 0, 0, t56 * t73, t58 * t73, -0.2e1 * t15 * t58 + 0.2e1 * (t20 + 0.2e1 * t44) * t77, 0.2e1 * t16 * t58 + 0.2e1 * (-t21 + 0.2e1 * t45) * t77, 0.2e1 * (-t15 * t54 - t16 * t53) * t56 + 0.2e1 * (-t20 * t54 - t21 * t53) * t76, 0.2e1 * pkin(6) ^ 2 * t71 + 0.2e1 * t20 * t15 + 0.2e1 * t21 * t16, 0.2e1 * t14 * t95 + 0.2e1 * t12 * t58 + 0.2e1 * (-t19 * t56 + t23 * t94) * qJD(2), 0.2e1 * (t12 * t54 - t53 * t9) * t56 + 0.2e1 * (-t18 * t53 + t19 * t54) * t76, -0.2e1 * t14 * t92 - 0.2e1 * t9 * t58 + 0.2e1 * (t18 * t56 - t23 * t91) * qJD(2), 0.2e1 * t19 * t12 + 0.2e1 * t23 * t14 + 0.2e1 * t18 * t9, 0.2e1 * t25 * t8, 0.2e1 * t8 * t24 - 0.2e1 * t25 * t7, -0.2e1 * t25 * t77 + 0.2e1 * t8 * t58, -0.2e1 * t24 * t77 - 0.2e1 * t7 * t58, -0.2e1 * t71, -0.2e1 * t10 * t24 + 0.2e1 * t17 * t7 + 0.2e1 * t2 * t58 - 0.2e1 * t65 * t77, 0.2e1 * t1 * t58 + 0.2e1 * t10 * t25 + 0.2e1 * t17 * t8 + 0.2e1 * t64 * t77; 0, 0, 0, 0, 0, t76, -t77, 0, -t49, t72, t40 + (t53 * t67 - t45) * qJD(2), t54 * t84 + (t54 * t67 + t44) * qJD(2), t63, -pkin(2) * t49 + (-t20 * t53 + t21 * t54) * qJD(3) + t63 * qJ(3), -t14 * t54 - t53 * t97 - t56 * t80 + t40, t12 * t53 + t9 * t54, -t14 * t53 + (t56 * t78 - t84 + t97) * t54, t14 * t34 + (qJ(3) * t9 + qJD(3) * t18) * t54 + (qJ(3) * t12 + qJD(3) * t19 - qJD(4) * t23) * t53, -t25 * t27 + t8 * t33, -t27 * t24 - t25 * t28 - t8 * t32 - t33 * t7, -t27 * t58 - t33 * t77, -t28 * t58 + t32 * t77, 0, t10 * t32 + t17 * t28 - t24 * t78 + t29 * t7 + t6 * t58 - t62 * t77, t10 * t33 - t17 * t27 + t25 * t78 + t29 * t8 + t5 * t58 + t61 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0.2e1 * t87, 0.2e1 * t54 * t78, t31, 0.2e1 * t80, -0.2e1 * t34 * t78 + 0.2e1 * t87, -0.2e1 * t33 * t27, 0.2e1 * t27 * t32 - 0.2e1 * t33 * t28, 0, 0, 0, 0.2e1 * t29 * t28 + 0.2e1 * t32 * t78, -0.2e1 * t29 * t27 + 0.2e1 * t33 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t41, 0, t49, t43, 0, -t41, t14, 0, 0, 0, 0, 0, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, 0, 0, 0, 0, 0, -t28, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t41, 0, t12, 0, 0, 0, 0, 0, -t55 * t81 - t57 * t77, t55 * t77 - t57 * t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, -t77, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t28, 0, t6, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, -t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t30;
