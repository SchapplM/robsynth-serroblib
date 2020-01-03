% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR7
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
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR7_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR7_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:36:22
% EndTime: 2019-12-31 19:36:25
% DurationCPUTime: 0.82s
% Computational Cost: add. (1180->113), mult. (2602->215), div. (0->0), fcn. (2365->6), ass. (0->80)
t40 = sin(qJ(5));
t37 = t40 ^ 2;
t42 = cos(qJ(5));
t38 = t42 ^ 2;
t65 = (t37 - t38) * qJD(5);
t89 = 2 * qJD(4);
t88 = pkin(3) + pkin(7);
t43 = cos(qJ(2));
t80 = cos(pkin(8));
t67 = t80 * t43;
t39 = sin(pkin(8));
t41 = sin(qJ(2));
t85 = t39 * t41;
t27 = -t67 + t85;
t82 = -qJ(3) - pkin(6);
t66 = qJD(2) * t82;
t51 = t43 * qJD(3) + t41 * t66;
t52 = -t41 * qJD(3) + t43 * t66;
t14 = t39 * t52 + t80 * t51;
t68 = t80 * t41;
t28 = t39 * t43 + t68;
t24 = t28 * qJD(2);
t8 = -t24 * pkin(4) + t14;
t87 = t27 * t8;
t86 = t37 * t24;
t23 = t38 * t24;
t84 = t40 * t24;
t83 = t42 * t24;
t79 = qJD(5) * t40;
t78 = qJD(5) * t42;
t77 = t41 * qJD(2);
t76 = t43 * qJD(2);
t75 = 0.2e1 * t27 * t24;
t25 = qJD(2) * t67 - t39 * t77;
t19 = 0.2e1 * t28 * t25;
t74 = -0.2e1 * pkin(1) * qJD(2);
t73 = t40 * t83;
t33 = t39 * pkin(2) + qJ(4);
t72 = t33 * t89;
t71 = pkin(2) * t77;
t70 = t41 * t76;
t69 = t40 * t78;
t36 = -t43 * pkin(2) - pkin(1);
t30 = t82 * t43;
t20 = -t39 * t30 - t82 * t68;
t26 = t27 ^ 2;
t64 = t26 * t69;
t60 = -t28 * qJ(4) + t36;
t11 = t88 * t27 + t60;
t59 = t28 * pkin(4) + t20;
t54 = t42 * t59;
t4 = -t40 * t11 + t54;
t5 = t42 * t11 + t40 * t59;
t63 = -t4 * t40 + t42 * t5;
t35 = -t80 * pkin(2) - pkin(3);
t13 = t39 * t51 - t80 * t52;
t21 = -t80 * t30 + t82 * t85;
t62 = t20 * t13 + t21 * t14;
t61 = -qJD(4) * t27 - t33 * t24;
t58 = t27 * t78 + t84;
t57 = t27 * t79 - t83;
t56 = t40 * t25 + t28 * t78;
t55 = -t42 * t25 + t28 * t79;
t53 = -0.2e1 * t28 * t24 - 0.2e1 * t25 * t27;
t32 = -pkin(7) + t35;
t50 = t8 + (t27 * t33 - t28 * t32) * qJD(5);
t49 = -t25 * qJ(4) - t28 * qJD(4) + t71;
t15 = -t27 * pkin(4) + t21;
t48 = -qJD(5) * t15 - t25 * t32 - t61;
t44 = t25 * pkin(4) + t13;
t45 = t88 * t24 + t49;
t2 = -qJD(5) * t54 + t11 * t79 - t40 * t44 - t42 * t45;
t3 = -qJD(5) * t5 - t40 * t45 + t42 * t44;
t1 = t63 * qJD(5) - t2 * t40 + t3 * t42;
t47 = -t2 * t42 - t3 * t40 + (-t4 * t42 - t40 * t5) * qJD(5);
t46 = 0.2e1 * t13 * t28 - 0.2e1 * t14 * t27 + 0.2e1 * t20 * t25 - 0.2e1 * t21 * t24;
t18 = t27 * pkin(3) + t60;
t10 = t24 * pkin(3) + t49;
t9 = -t27 * t65 + t73;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t70, 0.2e1 * (-t41 ^ 2 + t43 ^ 2) * qJD(2), 0, -0.2e1 * t70, 0, 0, t41 * t74, t43 * t74, 0, 0, t19, t53, 0, t75, 0, 0, 0.2e1 * t36 * t24 + 0.2e1 * t27 * t71, 0.2e1 * t36 * t25 + 0.2e1 * t28 * t71, t46, 0.2e1 * t36 * t71 + 0.2e1 * t62, 0, 0, 0, t19, t53, t75, t46, -0.2e1 * t10 * t27 - 0.2e1 * t18 * t24, -0.2e1 * t10 * t28 - 0.2e1 * t18 * t25, 0.2e1 * t18 * t10 + 0.2e1 * t62, 0.2e1 * t27 * t86 + 0.2e1 * t64, -0.2e1 * t26 * t65 + 0.4e1 * t27 * t73, 0.2e1 * t56 * t27 + 0.2e1 * t28 * t84, 0.2e1 * t27 * t23 - 0.2e1 * t64, -0.2e1 * t55 * t27 + 0.2e1 * t28 * t83, t19, 0.2e1 * t57 * t15 + 0.2e1 * t4 * t25 + 0.2e1 * t3 * t28 - 0.2e1 * t42 * t87, 0.2e1 * t58 * t15 + 0.2e1 * t2 * t28 - 0.2e1 * t5 * t25 + 0.2e1 * t40 * t87, 0.2e1 * t63 * t24 + 0.2e1 * t47 * t27, 0.2e1 * t15 * t8 - 0.2e1 * t5 * t2 + 0.2e1 * t4 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, 0, -t77, 0, -pkin(6) * t76, pkin(6) * t77, 0, 0, 0, 0, t25, 0, -t24, 0, -t13, -t14, (-t24 * t39 - t80 * t25) * pkin(2), (-t80 * t13 + t14 * t39) * pkin(2), 0, -t25, t24, 0, 0, 0, t35 * t25 + t61, t13, t14, t21 * qJD(4) + t13 * t35 + t14 * t33, t9, -0.4e1 * t27 * t69 + t23 - t86, -t55, -t9, -t56, 0, t50 * t40 - t48 * t42, t48 * t40 + t50 * t42, -t1, t15 * qJD(4) + t1 * t32 + t8 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t89, t72, -0.2e1 * t69, 0.2e1 * t65, 0, 0.2e1 * t69, 0, 0, 0.2e1 * qJD(4) * t40 + 0.2e1 * t33 * t78, 0.2e1 * qJD(4) * t42 - 0.2e1 * t33 * t79, 0, t72; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t25, 0, t71, 0, 0, 0, 0, 0, 0, 0, -t24, -t25, t10, 0, 0, 0, 0, 0, 0, -t56, t55, t23 + t86, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, 0, t13, 0, 0, 0, 0, 0, 0, -t55, -t56, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58, 0, -t57, t25, t3, t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, 0, -t78, 0, -t32 * t79, -t32 * t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t79, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t79, -t78, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
