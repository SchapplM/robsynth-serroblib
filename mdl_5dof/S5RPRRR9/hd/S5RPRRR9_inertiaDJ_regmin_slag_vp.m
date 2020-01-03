% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RPRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRRR9_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:59
% EndTime: 2019-12-31 19:08:01
% DurationCPUTime: 0.63s
% Computational Cost: add. (1273->104), mult. (2942->191), div. (0->0), fcn. (3027->8), ass. (0->79)
t54 = cos(qJ(5));
t48 = t54 ^ 2;
t51 = sin(qJ(5));
t83 = t51 ^ 2 - t48;
t73 = t83 * qJD(5);
t49 = sin(pkin(9));
t50 = cos(pkin(9));
t53 = sin(qJ(3));
t56 = cos(qJ(3));
t65 = t53 * t49 - t56 * t50;
t85 = pkin(6) + qJ(2);
t37 = t85 * t49;
t38 = t85 * t50;
t67 = -t56 * t37 - t53 * t38;
t94 = qJD(2) * t65 - qJD(3) * t67;
t41 = -t50 * pkin(2) - pkin(1);
t93 = 0.2e1 * t41;
t33 = t56 * t49 + t53 * t50;
t31 = t33 * qJD(3);
t92 = t31 * pkin(3);
t22 = -t33 * pkin(7) + t67;
t66 = t53 * t37 - t56 * t38;
t23 = -pkin(7) * t65 - t66;
t52 = sin(qJ(4));
t55 = cos(qJ(4));
t14 = t52 * t22 + t55 * t23;
t57 = -t31 * pkin(7) - t94;
t30 = t65 * qJD(3);
t59 = -qJD(2) * t33 + qJD(3) * t66;
t58 = -t30 * pkin(7) - t59;
t5 = qJD(4) * t14 + t52 * t57 + t55 * t58;
t3 = t5 * t51;
t13 = -t55 * t22 + t52 * t23;
t44 = qJD(5) * t54;
t91 = t13 * t44 + t3;
t68 = -t52 * t33 - t55 * t65;
t16 = qJD(4) * t68 - t55 * t30 - t52 * t31;
t26 = t55 * t33 - t52 * t65;
t90 = t26 * t16;
t89 = t26 * t54;
t17 = qJD(4) * t26 - t52 * t30 + t55 * t31;
t88 = t51 * t17;
t87 = t54 * t16;
t86 = t54 * t17;
t43 = -t55 * pkin(3) - pkin(4);
t82 = qJD(4) * t52;
t77 = pkin(3) * t82;
t84 = t43 * t44 + t51 * t77;
t81 = qJD(4) * t55;
t80 = qJD(5) * t51;
t79 = pkin(4) * t80;
t78 = pkin(4) * t44;
t76 = pkin(3) * t81;
t75 = t51 * t44;
t74 = -0.4e1 * t51 * t89;
t72 = 0.2e1 * (t49 ^ 2 + t50 ^ 2) * qJD(2);
t27 = pkin(3) * t65 + t41;
t15 = -pkin(4) * t68 - t26 * pkin(8) + t27;
t71 = t54 * t14 + t51 * t15;
t70 = t51 * t14 - t54 * t15;
t42 = t52 * pkin(3) + pkin(8);
t69 = -t26 * t43 - t42 * t68;
t64 = t43 * t80 - t54 * t77;
t63 = t51 * t16 + t26 * t44;
t62 = t26 * t80 - t87;
t10 = -t44 * t68 + t88;
t61 = -t68 * t80 - t86;
t60 = t16 * t43 - t17 * t42 + (t26 * t52 + t55 * t68) * qJD(4) * pkin(3);
t40 = 0.2e1 * t75;
t34 = -0.2e1 * t73;
t24 = t26 ^ 2;
t11 = t13 * t80;
t8 = -t26 * t73 + t51 * t87;
t7 = t17 * pkin(4) - t16 * pkin(8) + t92;
t6 = qJD(5) * t74 - t83 * t16;
t4 = -t22 * t81 + t23 * t82 + t52 * t58 - t55 * t57;
t2 = -t71 * qJD(5) + t51 * t4 + t54 * t7;
t1 = t70 * qJD(5) + t54 * t4 - t51 * t7;
t9 = [0, 0, 0, 0, 0, t72, qJ(2) * t72, -0.2e1 * t33 * t30, 0.2e1 * t30 * t65 - 0.2e1 * t33 * t31, 0, 0, 0, t31 * t93, -t30 * t93, 0.2e1 * t90, 0.2e1 * t16 * t68 - 0.2e1 * t26 * t17, 0, 0, 0, 0.2e1 * t27 * t17 - 0.2e1 * t68 * t92, 0.2e1 * t27 * t16 + 0.2e1 * t26 * t92, -0.2e1 * t24 * t75 + 0.2e1 * t48 * t90, t16 * t74 + 0.2e1 * t24 * t73, 0.2e1 * t26 * t86 + 0.2e1 * t62 * t68, -0.2e1 * t26 * t88 + 0.2e1 * t63 * t68, -0.2e1 * t68 * t17, 0.2e1 * t63 * t13 - 0.2e1 * t70 * t17 - 0.2e1 * t2 * t68 + 0.2e1 * t26 * t3, -0.2e1 * t1 * t68 - 0.2e1 * t62 * t13 - 0.2e1 * t71 * t17 + 0.2e1 * t5 * t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t30, 0, 0, 0, 0, 0, t17, t16, 0, 0, 0, 0, 0, -t61, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t31, 0, t59, t94, 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t61, 0, t11 + (-qJD(5) * t69 - t5) * t54 + t60 * t51, t54 * t60 + t69 * t80 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t77, -0.2e1 * t76, t40, t34, 0, 0, 0, 0.2e1 * t64, 0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, 0, -t5, t4, t8, t6, t10, -t61, 0, t11 + (-pkin(4) * t16 - pkin(8) * t17) * t51 + (-t5 + (-pkin(4) * t26 + pkin(8) * t68) * qJD(5)) * t54, pkin(4) * t62 + pkin(8) * t61 + t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, -t76, t40, t34, 0, 0, 0, t64 - t79, -t78 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, t34, 0, 0, 0, -0.2e1 * t79, -0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t63, t17, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t80, 0, -t42 * t44 - t51 * t76, t42 * t80 - t54 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t80, 0, -pkin(8) * t44, pkin(8) * t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
