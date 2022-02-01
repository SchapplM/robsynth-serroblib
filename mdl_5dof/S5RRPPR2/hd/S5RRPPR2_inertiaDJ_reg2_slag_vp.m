% Calculate inertial parameters regressor of joint inertia matrix time derivative for
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
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:51
% EndTime: 2022-01-20 10:05:54
% DurationCPUTime: 0.66s
% Computational Cost: add. (702->98), mult. (1592->174), div. (0->0), fcn. (1217->8), ass. (0->81)
t49 = sin(pkin(9));
t51 = cos(pkin(9));
t100 = -t51 * pkin(4) - t49 * pkin(7);
t99 = 2 * qJD(5);
t50 = sin(pkin(8));
t56 = cos(qJ(2));
t79 = pkin(1) * qJD(2);
t52 = cos(pkin(8));
t54 = sin(qJ(2));
t87 = t52 * t54;
t28 = (t50 * t56 + t87) * t79;
t55 = cos(qJ(5));
t66 = t56 * pkin(1) + pkin(2);
t81 = pkin(1) * t87 + t50 * t66;
t27 = qJ(4) + t81;
t53 = sin(qJ(5));
t61 = -t50 * t54 * pkin(1) + t52 * t66;
t59 = -pkin(3) - t61;
t58 = t59 + t100;
t88 = t51 * t55;
t6 = t27 * t88 + t53 * t58;
t78 = qJD(5) * t6;
t71 = t56 * t79;
t72 = t54 * t79;
t29 = -t50 * t72 + t52 * t71;
t25 = qJD(4) + t29;
t86 = t53 * t25;
t3 = t55 * t28 - t51 * t86 - t78;
t73 = t53 * qJD(4);
t44 = t50 * pkin(2) + qJ(4);
t65 = -t52 * pkin(2) - pkin(3);
t62 = -t65 - t100;
t13 = t44 * t88 - t53 * t62;
t77 = qJD(5) * t13;
t9 = -t51 * t73 - t77;
t98 = -t3 - t9;
t76 = qJD(5) * t53;
t41 = t51 * t76;
t57 = t55 * t58;
t2 = -qJD(5) * t57 - t25 * t88 + t27 * t41 - t53 * t28;
t97 = t2 * t53;
t60 = t55 * t62;
t74 = t51 * qJD(4);
t8 = qJD(5) * t60 + t44 * t41 - t55 * t74;
t94 = t8 * t53;
t47 = t49 ^ 2;
t21 = t47 * t25;
t93 = -t2 * t51 + t55 * t21;
t45 = t47 * qJD(4);
t92 = t55 * t45 - t8 * t51;
t91 = t28 * t49;
t90 = t28 * t51;
t48 = t51 ^ 2;
t22 = t48 * t25;
t89 = t51 * t53;
t85 = t44 * t21 + t27 * t45;
t75 = qJD(5) * t55;
t69 = t47 * t75;
t84 = t27 * t69 + t47 * t86;
t83 = t22 + t21;
t82 = t44 * t69 + t47 * t73;
t46 = t48 * qJD(4);
t80 = t46 + t45;
t70 = t47 * t76;
t68 = t49 * t76;
t67 = t49 * t75;
t64 = t49 * t51 * t99;
t63 = t53 * t69;
t42 = t51 * t75;
t35 = -0.2e1 * t63;
t34 = 0.2e1 * t63;
t33 = t55 * t64;
t32 = t53 * t64;
t31 = t44 * t45;
t26 = (t53 ^ 2 - t55 ^ 2) * t47 * t99;
t12 = -t44 * t89 - t60;
t11 = t12 * t68;
t10 = t27 * t21;
t5 = -t27 * t89 + t57;
t4 = t5 * t68;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t72, -0.2e1 * t71, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t28, -0.2e1 * t29, 0, -0.2e1 * t28 * t61 + 0.2e1 * t29 * t81, 0, 0, 0, 0, 0, 0, -0.2e1 * t90, 0.2e1 * t91, 0.2e1 * t83, 0.2e1 * t27 * t22 + 0.2e1 * t28 * t59 + 0.2e1 * t10, t35, t26, t32, t34, t33, 0, -0.2e1 * t3 * t51 + 0.2e1 * t84, -0.2e1 * t27 * t70 + 0.2e1 * t93, 0.2e1 * t4 + 0.2e1 * (t97 + (-t3 - t78) * t55) * t49, -0.2e1 * t6 * t2 + 0.2e1 * t5 * t3 + 0.2e1 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t72, -t71, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t29, 0, (-t28 * t52 + t29 * t50) * pkin(2), 0, 0, 0, 0, 0, 0, -t90, t91, t80 + t83, t28 * t65 + (t27 * qJD(4) + t25 * t44) * t48 + t85, t35, t26, t32, t34, t33, 0, t98 * t51 + t82 + t84, (-t27 - t44) * t70 + t92 + t93, t11 + t4 + ((t2 + t8) * t53 + ((-t13 - t6) * qJD(5) + t98) * t55) * t49, t3 * t12 - t2 * t13 + t5 * t9 - t6 * t8 + t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t80, 0.2e1 * t44 * t46 + 0.2e1 * t31, t35, t26, t32, t34, t33, 0, -0.2e1 * t9 * t51 + 0.2e1 * t82, -0.2e1 * t44 * t70 + 0.2e1 * t92, 0.2e1 * t11 + 0.2e1 * (t94 + (-t9 - t77) * t55) * t49, 0.2e1 * t12 * t9 - 0.2e1 * t13 * t8 + 0.2e1 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t2 * t55 - t51 * t25 - t3 * t53 + (-t5 * t55 - t53 * t6) * qJD(5)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (-t74 - t53 * t9 - t55 * t8 + (-t12 * t55 - t13 * t53) * qJD(5)) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, t41, t42, 0, -t97 + t3 * t55 + (-t5 * t53 + t55 * t6) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t42, 0, -t94 + t9 * t55 + (-t12 * t53 + t13 * t55) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, -t67, 0, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, 0, -t67, 0, t9, t8, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t67, t68, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
