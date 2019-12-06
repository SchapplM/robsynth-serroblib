% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRP7_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:56:26
% EndTime: 2019-12-05 16:56:29
% DurationCPUTime: 0.70s
% Computational Cost: add. (483->131), mult. (1410->256), div. (0->0), fcn. (1189->8), ass. (0->85)
t38 = sin(qJ(3));
t92 = -0.4e1 * t38;
t41 = cos(qJ(3));
t51 = -t41 * pkin(3) - t38 * pkin(8);
t25 = -pkin(2) + t51;
t40 = cos(qJ(4));
t87 = t40 * t41;
t29 = pkin(7) * t87;
t37 = sin(qJ(4));
t83 = t37 * t25 + t29;
t33 = t40 ^ 2;
t82 = t37 ^ 2 - t33;
t55 = t82 * qJD(4);
t91 = pkin(7) * t37;
t35 = sin(pkin(5));
t39 = sin(qJ(2));
t90 = t35 * t39;
t42 = cos(qJ(2));
t89 = t35 * t42;
t88 = t38 * t40;
t86 = -qJ(5) - pkin(8);
t50 = pkin(3) * t38 - pkin(8) * t41;
t23 = t50 * qJD(3);
t74 = qJD(4) * t40;
t85 = -t37 * t23 - t25 * t74;
t72 = t38 * qJD(3);
t60 = t37 * t72;
t84 = pkin(7) * t60 + t40 * t23;
t32 = t38 ^ 2;
t81 = -t41 ^ 2 + t32;
t80 = qJ(5) * t38;
t79 = t40 * qJ(5);
t78 = qJD(2) * t39;
t36 = cos(pkin(5));
t16 = -t36 * t41 + t38 * t90;
t77 = qJD(3) * t16;
t76 = qJD(3) * t40;
t75 = qJD(4) * t37;
t73 = qJD(4) * t41;
t71 = t40 * qJD(5);
t70 = t41 * qJD(3);
t69 = -0.2e1 * pkin(2) * qJD(3);
t68 = -0.2e1 * pkin(3) * qJD(4);
t67 = pkin(4) * t75;
t66 = pkin(7) * t70;
t65 = t37 * t73;
t64 = t40 * t73;
t63 = t16 * t75;
t62 = t35 * t78;
t61 = qJD(2) * t89;
t59 = t37 * t74;
t58 = t38 * t70;
t57 = t40 * t70;
t56 = qJD(4) * t86;
t54 = t81 * qJD(3);
t53 = 0.2e1 * t58;
t52 = t37 * t57;
t17 = t36 * t38 + t41 * t90;
t10 = -t17 * t37 - t40 * t89;
t48 = -t17 * t40 + t37 * t89;
t49 = -t10 * t40 + t37 * t48;
t9 = t17 * qJD(3) + t38 * t61;
t47 = t16 * t74 + t9 * t37;
t46 = -t9 * t40 + t63;
t45 = -t38 * t75 + t57;
t44 = t40 * t72 + t65;
t43 = t37 * t70 + t38 * t74;
t30 = -t40 * pkin(4) - pkin(3);
t27 = t86 * t40;
t26 = t86 * t37;
t24 = (pkin(4) * t37 + pkin(7)) * t38;
t22 = t40 * t25;
t15 = -t37 * qJD(5) + t40 * t56;
t14 = t37 * t56 + t71;
t13 = t43 * pkin(4) + t66;
t12 = -t37 * t80 + t83;
t8 = t41 * t61 - t77;
t7 = -t38 * t79 + t22 + (-pkin(4) - t91) * t41;
t6 = -t83 * qJD(4) + t84;
t5 = t44 * pkin(7) + t85;
t4 = (-pkin(7) * qJD(3) - qJ(5) * qJD(4)) * t88 + (-qJD(5) * t38 + (-pkin(7) * qJD(4) - qJ(5) * qJD(3)) * t41) * t37 - t85;
t3 = qJD(4) * t10 + t37 * t62 + t8 * t40;
t2 = t48 * qJD(4) - t8 * t37 + t40 * t62;
t1 = -t38 * t71 + (pkin(4) * t38 - t41 * t79) * qJD(3) + (-t29 + (-t25 + t80) * t37) * qJD(4) + t84;
t11 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t10 * t2 + 0.2e1 * t16 * t9 - 0.2e1 * t3 * t48; 0, 0, -t62, -t61, 0, 0, 0, 0, 0, (-t41 * t78 - t42 * t72) * t35, (t38 * t78 - t42 * t70) * t35, 0, 0, 0, 0, 0, (t37 * t77 - t2) * t41 + (qJD(3) * t10 + t47) * t38, (t16 * t76 + t3) * t41 + (qJD(3) * t48 - t46) * t38, t49 * t70 + (-t2 * t40 - t3 * t37 + (t10 * t37 + t40 * t48) * qJD(4)) * t38, t10 * t1 + t3 * t12 + t16 * t13 + t2 * t7 + t9 * t24 - t4 * t48; 0, 0, 0, 0, t53, -0.2e1 * t54, 0, 0, 0, t38 * t69, t41 * t69, -0.2e1 * t32 * t59 + 0.2e1 * t33 * t58, 0.2e1 * t32 * t55 + t52 * t92, 0.2e1 * t38 * t65 + 0.2e1 * t81 * t76, -0.2e1 * t37 * t54 + 0.2e1 * t38 * t64, -0.2e1 * t58, 0.2e1 * t22 * t72 - 0.2e1 * t6 * t41 + 0.2e1 * (t32 * t74 + t37 * t58) * pkin(7), -0.2e1 * t5 * t41 - 0.2e1 * t83 * t72 + 0.2e1 * (-t32 * t75 + t40 * t53) * pkin(7), 0.2e1 * (-t12 * t37 - t40 * t7) * t70 + 0.2e1 * (-t1 * t40 - t37 * t4 + (-t12 * t40 + t37 * t7) * qJD(4)) * t38, 0.2e1 * t7 * t1 + 0.2e1 * t12 * t4 + 0.2e1 * t24 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t8, 0, 0, 0, 0, 0, t46, t47, t49 * qJD(4) - t2 * t37 + t3 * t40, pkin(4) * t63 + t10 * t15 - t14 * t48 + t2 * t26 - t3 * t27 + t9 * t30; 0, 0, 0, 0, 0, 0, t70, -t72, 0, -t66, pkin(7) * t72, -t38 * t55 + t52, t59 * t92 - t82 * t70, t60 - t64, t44, 0, (pkin(8) * t87 + (-pkin(3) * t40 + t91) * t38) * qJD(4) + (t51 * t37 - t29) * qJD(3), (pkin(7) * t88 + t50 * t37) * qJD(4) + (t51 * t40 + t41 * t91) * qJD(3), (-t26 * t70 - t15 * t38 + t4 + (t27 * t38 - t7) * qJD(4)) * t40 + (t27 * t70 - t14 * t38 - t1 + (t26 * t38 - t12) * qJD(4)) * t37, t1 * t26 + t12 * t14 + t13 * t30 + t7 * t15 + t24 * t67 - t4 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t59, -0.2e1 * t55, 0, 0, 0, t37 * t68, t40 * t68, 0.2e1 * t14 * t40 - 0.2e1 * t15 * t37 + 0.2e1 * (-t26 * t40 + t27 * t37) * qJD(4), -0.2e1 * t27 * t14 + 0.2e1 * t26 * t15 + 0.2e1 * t30 * t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t3, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t43, t72, t6, t5, -t45 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, -t75, 0, -pkin(8) * t74, pkin(8) * t75, -pkin(4) * t74, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t67; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t11;
