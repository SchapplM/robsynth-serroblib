% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RRPPR5
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
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRPPR5_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR5_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:29:58
% EndTime: 2019-12-31 19:30:00
% DurationCPUTime: 0.65s
% Computational Cost: add. (975->104), mult. (2154->193), div. (0->0), fcn. (1928->6), ass. (0->62)
t56 = cos(qJ(2));
t55 = sin(qJ(2));
t81 = sin(pkin(8));
t69 = t81 * t55;
t82 = cos(pkin(8));
t41 = -t82 * t56 + t69;
t85 = 2 * qJD(4);
t84 = cos(qJ(5));
t83 = -qJ(3) - pkin(6);
t47 = t83 * t56;
t27 = -t82 * t47 + t83 * t69;
t54 = sin(qJ(5));
t80 = qJD(5) * t54;
t79 = t55 * qJD(2);
t78 = t56 * qJD(2);
t71 = t82 * t55;
t42 = t81 * t56 + t71;
t38 = t42 * qJD(2);
t77 = 0.2e1 * t41 * t38;
t76 = -0.2e1 * pkin(1) * qJD(2);
t75 = pkin(2) * t79;
t74 = t55 * t78;
t52 = -t56 * pkin(2) - pkin(1);
t68 = qJD(2) * t83;
t62 = t56 * qJD(3) + t55 * t68;
t63 = -t55 * qJD(3) + t56 * t68;
t17 = t81 * t62 - t82 * t63;
t18 = t82 * t62 + t81 * t63;
t26 = -t81 * t47 - t83 * t71;
t73 = t26 * t17 + t27 * t18;
t72 = qJD(5) * t84;
t51 = -t82 * pkin(2) - pkin(3);
t39 = t41 * qJD(2);
t67 = t42 * t38 - t39 * t41;
t66 = -pkin(4) + t51;
t65 = t42 * qJ(4) - t52;
t22 = t54 * t41 + t84 * t42;
t12 = t38 * pkin(3) + t39 * qJ(4) - t42 * qJD(4) + t75;
t64 = -t42 * pkin(7) + t26;
t61 = t84 * t66;
t60 = 0.2e1 * t17 * t42 - 0.2e1 * t18 * t41 - 0.2e1 * t26 * t39 - 0.2e1 * t27 * t38;
t59 = t84 * t64;
t49 = t81 * pkin(2) + qJ(4);
t29 = t84 * t49 + t54 * t66;
t19 = t41 * pkin(7) + t27;
t4 = t84 * t19 + t54 * t64;
t58 = t39 * pkin(7) + t17;
t57 = t38 * pkin(7) + t18;
t28 = -t54 * t49 + t61;
t25 = t54 * qJD(4) + t29 * qJD(5);
t24 = -t84 * qJD(4) - qJD(5) * t61 + t49 * t80;
t23 = -0.2e1 * t42 * t39;
t21 = -t84 * t41 + t54 * t42;
t20 = t41 * pkin(3) - t65;
t14 = (-pkin(3) - pkin(4)) * t41 + t65;
t7 = t38 * pkin(4) + t12;
t6 = qJD(5) * t22 - t84 * t38 - t54 * t39;
t5 = -t54 * t38 + t84 * t39 - t41 * t72 + t42 * t80;
t3 = -t54 * t19 + t59;
t2 = t4 * qJD(5) + t54 * t57 - t84 * t58;
t1 = -qJD(5) * t59 + t19 * t80 - t54 * t58 - t84 * t57;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t74, 0.2e1 * (-t55 ^ 2 + t56 ^ 2) * qJD(2), 0, -0.2e1 * t74, 0, 0, t55 * t76, t56 * t76, 0, 0, t23, -0.2e1 * t67, 0, t77, 0, 0, 0.2e1 * t52 * t38 + 0.2e1 * t41 * t75, -0.2e1 * t52 * t39 + 0.2e1 * t42 * t75, t60, 0.2e1 * t52 * t75 + 0.2e1 * t73, t23, 0, 0.2e1 * t67, 0, 0, t77, 0.2e1 * t12 * t41 + 0.2e1 * t20 * t38, t60, -0.2e1 * t12 * t42 + 0.2e1 * t20 * t39, 0.2e1 * t20 * t12 + 0.2e1 * t73, -0.2e1 * t22 * t5, 0.2e1 * t5 * t21 - 0.2e1 * t22 * t6, 0, 0.2e1 * t21 * t6, 0, 0, 0.2e1 * t14 * t6 - 0.2e1 * t7 * t21, -0.2e1 * t14 * t5 - 0.2e1 * t7 * t22, 0.2e1 * t1 * t21 + 0.2e1 * t2 * t22 + 0.2e1 * t3 * t5 - 0.2e1 * t4 * t6, -0.2e1 * t4 * t1 - 0.2e1 * t14 * t7 - 0.2e1 * t3 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, 0, -t79, 0, -pkin(6) * t78, pkin(6) * t79, 0, 0, 0, 0, -t39, 0, -t38, 0, -t17, -t18, (-t81 * t38 + t82 * t39) * pkin(2), (-t82 * t17 + t81 * t18) * pkin(2), 0, -t39, 0, 0, t38, 0, -t17, -qJD(4) * t41 - t49 * t38 - t51 * t39, t18, t27 * qJD(4) + t17 * t51 + t18 * t49, 0, 0, t5, 0, t6, 0, t2, -t1, t24 * t21 + t25 * t22 + t28 * t5 - t29 * t6, -t1 * t29 - t2 * t28 - t4 * t24 - t3 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85, t49 * t85, 0, 0, 0, 0, 0, 0, 0.2e1 * t25, -0.2e1 * t24, 0, -0.2e1 * t29 * t24 - 0.2e1 * t28 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, t75, 0, 0, 0, 0, 0, 0, t38, 0, t39, t12, 0, 0, 0, 0, 0, 0, -t6, t5, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, 0, t17, 0, 0, 0, 0, 0, 0, 0, 0, t84 * t5 - t54 * t6 + (-t84 * t21 + t22 * t54) * qJD(5), -t2 * t84 - t1 * t54 + (-t3 * t54 + t84 * t4) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t72, 0, -t25 * t84 - t24 * t54 + (-t28 * t54 + t84 * t29) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, 0, -t6, 0, -t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t24, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t80, -t72, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t8;
