% Calculate minimal parameter regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 18:09
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRP6_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaDJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 18:08:28
% EndTime: 2021-01-15 18:08:32
% DurationCPUTime: 0.69s
% Computational Cost: add. (564->129), mult. (1353->234), div. (0->0), fcn. (984->6), ass. (0->87)
t40 = sin(qJ(3));
t97 = -0.4e1 * t40;
t41 = cos(qJ(4));
t33 = qJD(4) * t41;
t39 = sin(qJ(4));
t42 = cos(qJ(3));
t72 = t42 * qJD(3);
t16 = t40 * t33 + t39 * t72;
t96 = t16 * pkin(4);
t28 = sin(pkin(8)) * pkin(1) + pkin(6);
t95 = t40 * qJD(5) + (qJ(5) * qJD(3) + qJD(4) * t28) * t42;
t94 = 0.2e1 * qJD(4);
t93 = t41 * pkin(4);
t58 = t28 * t72;
t7 = t58 + t96;
t92 = t7 * t39;
t91 = t7 * t41;
t29 = -cos(pkin(8)) * pkin(1) - pkin(2);
t53 = -t42 * pkin(3) - t40 * pkin(7);
t46 = t29 + t53;
t45 = qJD(4) * t46;
t52 = pkin(3) * t40 - pkin(7) * t42;
t47 = t52 * qJD(3);
t90 = -t39 * t47 - t41 * t45;
t83 = t41 * t42;
t21 = t28 * t83;
t9 = t39 * t46;
t89 = t21 + t9;
t82 = -qJ(5) - pkin(7);
t23 = t82 * t39;
t88 = t23 * t40;
t24 = t82 * t41;
t87 = t24 * t40;
t86 = t28 * t39;
t85 = t39 * t42;
t84 = t40 * t41;
t31 = t40 * qJD(3);
t60 = t28 * t31;
t81 = t39 * t60 + t41 * t47;
t34 = t39 ^ 2;
t36 = t41 ^ 2;
t80 = t34 - t36;
t35 = t40 ^ 2;
t79 = -t42 ^ 2 + t35;
t78 = qJ(5) * t40;
t13 = (pkin(4) * t39 + t28) * t40;
t77 = qJD(3) * t13;
t76 = qJD(3) * t41;
t75 = qJD(4) * t39;
t74 = qJD(4) * t42;
t71 = -0.2e1 * pkin(3) * qJD(4);
t70 = 0.2e1 * qJD(3) * t29;
t69 = pkin(4) * t31;
t68 = pkin(4) * t75;
t67 = t40 * t75;
t66 = t39 * t74;
t65 = t41 * t74;
t64 = t13 * t75;
t63 = t34 * t72;
t62 = t39 * t33;
t61 = t40 * t72;
t59 = t41 * t72;
t30 = -pkin(3) - t93;
t57 = -t30 + t93;
t56 = t80 * qJD(4);
t55 = t79 * qJD(3);
t54 = t39 * t59;
t10 = t41 * t46;
t5 = -t41 * t78 + t10 + (-pkin(4) - t86) * t42;
t6 = -t39 * t78 + t89;
t51 = -t39 * t6 - t41 * t5;
t50 = t39 * t5 - t41 * t6;
t49 = pkin(4) * t34 + t30 * t41;
t14 = -t59 + t67;
t15 = t41 * t31 + t66;
t11 = -t41 * qJD(5) - t82 * t75;
t12 = -t39 * qJD(5) + t82 * t33;
t44 = -t11 * t41 - t12 * t39 + (-t23 * t41 + t24 * t39) * qJD(4);
t43 = qJ(5) * t67 - t39 * t45 - t95 * t41 + t81;
t26 = t36 * t72;
t22 = t36 * t61;
t17 = t39 * t31 - t65;
t4 = -t89 * qJD(4) + t81;
t3 = t15 * t28 + t90;
t2 = (qJ(5) * qJD(4) + qJD(3) * t28) * t84 + t95 * t39 + t90;
t1 = t43 + t69;
t8 = [0, 0, 0, 0, 0.2e1 * t61, -0.2e1 * t55, 0, 0, 0, t40 * t70, t42 * t70, -0.2e1 * t35 * t62 + 0.2e1 * t22, t80 * t35 * t94 + t54 * t97, 0.2e1 * t40 * t66 + 0.2e1 * t79 * t76, -0.2e1 * t39 * t55 + 0.2e1 * t40 * t65, -0.2e1 * t61, 0.2e1 * t10 * t31 - 0.2e1 * t4 * t42 + 0.2e1 * (t35 * t33 + t39 * t61) * t28, -0.2e1 * t35 * t28 * t75 - 0.2e1 * t3 * t42 + 0.2e1 * (t21 - t9) * t31, 0.2e1 * (t39 * t77 - t1) * t42 + 0.2e1 * (qJD(3) * t5 + t13 * t33 + t92) * t40, 0.2e1 * (t13 * t76 - t2) * t42 + 0.2e1 * (-qJD(3) * t6 - t64 + t91) * t40, 0.2e1 * t51 * t72 + 0.2e1 * (t50 * qJD(4) - t1 * t41 + t2 * t39) * t40, 0.2e1 * t5 * t1 + 0.2e1 * t13 * t7 - 0.2e1 * t6 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t50 * qJD(3) - t7) * t42 + (t51 * qJD(4) - t1 * t39 - t2 * t41 + t77) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t22 + 0.2e1 * (t34 - 0.1e1) * t61; 0, 0, 0, 0, 0, 0, t72, -t31, 0, -t58, t60, -t40 * t56 + t54, t62 * t97 + t26 - t63, t17, t15, 0, (pkin(7) * t83 + (-pkin(3) * t41 + t86) * t40) * qJD(4) + (t53 * t39 - t21) * qJD(3), (t28 * t84 + t52 * t39) * qJD(4) + (t28 * t85 + t53 * t41) * qJD(3), -t12 * t42 - t91 + (t30 * t85 + t88) * qJD(3) + (t13 * t39 + t49 * t40) * qJD(4), -t11 * t42 + t92 + (t30 * t83 + t87) * qJD(3) + (t57 * t40 * t39 + t13 * t41) * qJD(4), (-t23 * t72 - t12 * t40 - t2 + (-t5 + t87) * qJD(4)) * t41 + (t24 * t72 + t11 * t40 - t1 + (-t6 + t88) * qJD(4)) * t39, pkin(4) * t64 + t1 * t23 - t6 * t11 + t5 * t12 + t2 * t24 + t7 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t72, 0, 0, 0, 0, 0, -t15, t17, -t15, t17, t26 + t63, (-t68 + (-t23 * t39 - t24 * t41) * qJD(3)) * t42 + (qJD(3) * t30 + t44) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t62, -0.2e1 * t56, 0, 0, 0, t39 * t71, t41 * t71, -0.2e1 * t57 * t75, t49 * t94, 0.2e1 * t44, 0.2e1 * t24 * t11 + 0.2e1 * t23 * t12 + 0.2e1 * t30 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t16, t31, t4, t3, t43 + 0.2e1 * t69, t2, t14 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t14, -t16, t14, 0, -t96; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t75, 0, -pkin(7) * t33, pkin(7) * t75, t12, t11, -pkin(4) * t33, t12 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t33, 0, t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
