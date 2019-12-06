% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:40
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR2_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR2_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR2_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:39:56
% EndTime: 2019-12-05 17:40:00
% DurationCPUTime: 0.63s
% Computational Cost: add. (1001->77), mult. (1954->145), div. (0->0), fcn. (1926->6), ass. (0->57)
t44 = sin(pkin(8));
t45 = cos(pkin(8));
t48 = sin(qJ(4));
t49 = cos(qJ(4));
t34 = t49 * t44 + t48 * t45;
t30 = t34 * qJD(4);
t47 = sin(qJ(5));
t69 = t48 * t44;
t35 = t49 * t45 - t69;
t72 = cos(qJ(5));
t54 = t47 * t34 - t72 * t35;
t66 = qJD(4) * t49;
t67 = qJD(4) * t48;
t58 = t44 * t67 - t45 * t66;
t5 = t54 * qJD(5) + t47 * t30 + t72 * t58;
t55 = t72 * t34 + t47 * t35;
t78 = t5 * t55;
t56 = t34 * t58;
t71 = t35 * t30;
t88 = 0.2e1 * t56 + 0.2e1 * t71;
t46 = -pkin(1) - qJ(3);
t73 = -pkin(6) + t46;
t36 = t73 * t44;
t62 = t73 * t45;
t20 = t49 * t36 + t48 * t62;
t14 = -t34 * pkin(7) + t20;
t13 = qJD(3) * t69 - t36 * t66 + (-t49 * qJD(3) - t73 * t67) * t45;
t50 = t30 * pkin(7) + t13;
t33 = t49 * t62;
t12 = t34 * qJD(3) - qJD(4) * t33 + t36 * t67;
t51 = -t58 * pkin(7) + t12;
t19 = -t48 * t36 + t33;
t57 = t35 * pkin(7) - t19;
t52 = t72 * t57;
t65 = qJD(5) * t47;
t1 = qJD(5) * t52 + t14 * t65 - t47 * t50 + t72 * t51;
t4 = t72 * t14 - t47 * t57;
t2 = -t4 * qJD(5) + t47 * t51 + t72 * t50;
t87 = t1 * t55 + t2 * t54 + t4 * t5;
t86 = -(t47 * t54 + t55 * t72) * qJD(5) + t47 * t5;
t61 = qJD(5) * t72;
t75 = -t72 * t30 + t47 * t58;
t6 = t34 * t61 + t35 * t65 - t75;
t85 = t54 * t6;
t7 = -t55 * qJD(5) + t75;
t82 = -t54 * t7 - t78;
t37 = (t44 ^ 2 + t45 ^ 2) * qJD(3);
t76 = t12 * t34 - t13 * t35 + t19 * t30 + t20 * t58;
t74 = 2 * qJD(2);
t39 = t44 * pkin(3) + qJ(2);
t64 = qJ(2) * qJD(2);
t63 = pkin(4) * t65;
t60 = pkin(4) * t61;
t22 = -t58 * pkin(4) + qJD(2);
t21 = t34 * pkin(4) + t39;
t3 = -t47 * t14 - t52;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t74, 0.2e1 * t64, 0, 0, 0, 0, 0, 0, t44 * t74, t45 * t74, 0.2e1 * t37, -0.2e1 * t46 * t37 + 0.2e1 * t64, -0.2e1 * t71, 0.2e1 * t30 * t34 + 0.2e1 * t35 * t58, 0, -0.2e1 * t56, 0, 0, 0.2e1 * qJD(2) * t34 - 0.2e1 * t39 * t58, 0.2e1 * qJD(2) * t35 - 0.2e1 * t39 * t30, 0.2e1 * t76, 0.2e1 * t39 * qJD(2) - 0.2e1 * t20 * t12 + 0.2e1 * t19 * t13, 0.2e1 * t85, -0.2e1 * t5 * t54 + 0.2e1 * t55 * t6, 0, -0.2e1 * t78, 0, 0, -0.2e1 * t21 * t5 + 0.2e1 * t22 * t55, -0.2e1 * t21 * t6 - 0.2e1 * t22 * t54, 0.2e1 * t3 * t6 + 0.2e1 * t87, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t21 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, 0, 0, 0, 0, 0, 0, 0, 0, t88, -t76, 0, 0, 0, 0, 0, 0, 0, 0, t78 - t82 - t85, t3 * t7 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t88, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, -t58, -t30, 0, qJD(2), 0, 0, 0, 0, 0, 0, -t5, -t6, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, 0, t58, 0, t13, t12, 0, 0, 0, 0, -t6, 0, t5, 0, t2, t1, (t72 * t6 + t86) * pkin(4), (t72 * t2 - t1 * t47 + (-t3 * t47 + t72 * t4) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t58, 0, 0, 0, 0, 0, 0, 0, 0, t7, t5, 0, (t72 * t7 - t86) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t63, -0.2e1 * t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, t5, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, -t60, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
