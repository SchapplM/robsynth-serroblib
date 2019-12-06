% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:55
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR5_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:54:41
% EndTime: 2019-12-05 15:54:45
% DurationCPUTime: 0.68s
% Computational Cost: add. (944->103), mult. (2415->202), div. (0->0), fcn. (2392->8), ass. (0->64)
t50 = sin(pkin(9));
t51 = cos(pkin(9));
t55 = cos(qJ(4));
t75 = qJD(4) * t55;
t53 = sin(qJ(4));
t76 = qJD(4) * t53;
t35 = t50 * t75 + t51 * t76;
t87 = 0.2e1 * t35;
t82 = t53 * t51;
t40 = t55 * t50 + t82;
t54 = sin(qJ(2));
t30 = t40 * t54;
t63 = t53 * t50 - t55 * t51;
t79 = pkin(6) + qJ(3);
t69 = t79 * t50;
t38 = t55 * t69;
t42 = t79 * t51;
t73 = t55 * qJD(3);
t12 = (t50 * qJD(3) + qJD(4) * t42) * t53 + qJD(4) * t38 - t51 * t73;
t34 = t63 * qJD(4);
t86 = -0.2e1 * t34;
t85 = t35 * pkin(4);
t84 = cos(qJ(5));
t56 = cos(qJ(2));
t80 = t56 * t35;
t24 = t55 * t42 - t53 * t69;
t77 = t50 ^ 2 + t51 ^ 2;
t52 = sin(qJ(5));
t74 = qJD(5) * t52;
t47 = t54 * qJD(2);
t72 = t56 * qJD(2);
t71 = pkin(4) * t74;
t70 = t54 * t72;
t46 = -t51 * pkin(3) - pkin(2);
t68 = t77 * t56;
t23 = -t53 * t42 - t38;
t67 = qJD(5) * t84;
t66 = t77 * qJD(3);
t65 = pkin(4) * t67;
t64 = 0.2e1 * t66;
t61 = t40 * pkin(7) - t23;
t31 = t63 * t54;
t17 = -t52 * t30 - t84 * t31;
t60 = t84 * t63;
t59 = t84 * t61;
t22 = t84 * t40 - t52 * t63;
t18 = -t63 * pkin(7) + t24;
t6 = t84 * t18 - t52 * t61;
t58 = t35 * pkin(7) + t12;
t13 = -qJD(3) * t82 - t42 * t75 + (t79 * t76 - t73) * t50;
t57 = t34 * pkin(7) + t13;
t26 = t63 * pkin(4) + t46;
t21 = t52 * t40 + t60;
t20 = t54 * t34 - t40 * t72;
t19 = qJD(4) * t30 + t63 * t72;
t16 = -t84 * t30 + t52 * t31;
t8 = qJD(5) * t22 - t52 * t34 + t84 * t35;
t7 = qJD(5) * t60 + t84 * t34 + t52 * t35 + t40 * t74;
t5 = -t52 * t18 - t59;
t4 = -t17 * qJD(5) + t52 * t19 + t84 * t20;
t3 = t84 * t19 - t52 * t20 + t30 * t67 - t31 * t74;
t2 = -t6 * qJD(5) + t52 * t58 + t84 * t57;
t1 = qJD(5) * t59 + t18 * t74 - t52 * t57 + t84 * t58;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * (-0.1e1 + t77) * t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t31 * t19 - 0.2e1 * t30 * t20 - 0.2e1 * t70, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t16 * t4 - 0.2e1 * t17 * t3 - 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, -t72, 0, 0, 0, 0, 0, 0, 0, 0, -t51 * t47, t50 * t47, qJD(2) * t68, t54 * t66 + (-pkin(2) * t54 + qJ(3) * t68) * qJD(2), 0, 0, 0, 0, 0, 0, t63 * t47 - t80, t56 * t34 + t40 * t47, t19 * t63 - t20 * t40 - t30 * t34 + t31 * t35, t31 * t12 - t30 * t13 - t19 * t24 + t20 * t23 + t46 * t47, 0, 0, 0, 0, 0, 0, t21 * t47 - t56 * t8, t22 * t47 + t56 * t7, t16 * t7 - t17 * t8 + t3 * t21 - t4 * t22, -pkin(4) * t80 - t17 * t1 + t16 * t2 + t26 * t47 - t3 * t6 + t4 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, qJ(3) * t64, t40 * t86, 0.2e1 * t34 * t63 - 0.2e1 * t40 * t35, 0, t63 * t87, 0, 0, t46 * t87, t46 * t86, 0.2e1 * t12 * t63 - 0.2e1 * t13 * t40 + 0.2e1 * t23 * t34 - 0.2e1 * t24 * t35, -0.2e1 * t24 * t12 + 0.2e1 * t23 * t13, -0.2e1 * t22 * t7, 0.2e1 * t7 * t21 - 0.2e1 * t22 * t8, 0, 0.2e1 * t21 * t8, 0, 0, 0.2e1 * t21 * t85 + 0.2e1 * t26 * t8, 0.2e1 * t22 * t85 - 0.2e1 * t26 * t7, 0.2e1 * t1 * t21 - 0.2e1 * t2 * t22 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t26 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47, 0, 0, 0, 0, 0, 0, 0, 0, 0, t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t34, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t19, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, (t84 * t4 - t3 * t52 + (-t16 * t52 + t84 * t17) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, 0, -t35, 0, t13, t12, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, (t84 * t7 - t52 * t8 + (-t84 * t21 + t22 * t52) * qJD(5)) * pkin(4), (t84 * t2 - t1 * t52 + (-t5 * t52 + t84 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t71, -0.2e1 * t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t65, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t9;
