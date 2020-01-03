% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPRPR4_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR4_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR4_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:39:17
% EndTime: 2020-01-03 11:39:22
% DurationCPUTime: 0.65s
% Computational Cost: add. (1070->92), mult. (2275->179), div. (0->0), fcn. (2148->8), ass. (0->59)
t48 = sin(qJ(5));
t73 = cos(pkin(9));
t61 = t73 * pkin(3) + pkin(4);
t72 = sin(pkin(9));
t66 = t72 * pkin(3);
t77 = cos(qJ(5));
t32 = -t48 * t66 + t77 * t61;
t80 = 2 * qJD(3);
t49 = sin(qJ(3));
t50 = cos(qJ(3));
t38 = t73 * t49 + t72 * t50;
t56 = t72 * t49 - t73 * t50;
t54 = t77 * t56;
t17 = t48 * t38 + t54;
t18 = t77 * t38 - t48 * t56;
t63 = qJD(3) * t72;
t64 = qJD(3) * t73;
t35 = -t49 * t63 + t50 * t64;
t75 = t49 * t64 + t50 * t63;
t8 = t18 * qJD(5) + t48 * t35 + t77 * t75;
t79 = t17 * t8;
t71 = qJD(5) * t48;
t7 = qJD(5) * t54 - t77 * t35 + t38 * t71 + t48 * t75;
t78 = t18 * t7;
t76 = t38 * t35;
t44 = sin(pkin(8)) * pkin(1) + pkin(6);
t74 = qJ(4) + t44;
t36 = t74 * t49;
t37 = t74 * t50;
t16 = -t72 * t36 + t73 * t37;
t70 = t49 * qJD(3);
t69 = t50 * qJD(3);
t45 = -cos(pkin(8)) * pkin(1) - pkin(2);
t68 = t45 * t80;
t46 = pkin(3) * t70;
t67 = t49 * t69;
t65 = qJD(3) * t74;
t15 = -t73 * t36 - t72 * t37;
t40 = -t50 * pkin(3) + t45;
t60 = t38 * pkin(7) - t15;
t58 = -t49 * qJD(4) - t50 * t65;
t57 = t50 * qJD(4) - t49 * t65;
t55 = t77 * t60;
t53 = t56 * t75;
t33 = t48 * t61 + t77 * t66;
t14 = -t56 * pkin(7) + t16;
t4 = t77 * t14 - t48 * t60;
t13 = t73 * t57 + t72 * t58;
t12 = -t72 * t57 + t73 * t58;
t52 = -t35 * pkin(7) + t12;
t51 = -t75 * pkin(7) + t13;
t24 = t33 * qJD(5);
t23 = t32 * qJD(5);
t22 = t75 * pkin(4) + t46;
t21 = t56 * pkin(4) + t40;
t3 = -t48 * t14 - t55;
t2 = -t4 * qJD(5) - t48 * t51 + t77 * t52;
t1 = qJD(5) * t55 + t14 * t71 - t48 * t52 - t77 * t51;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t67, (-t49 ^ 2 + t50 ^ 2) * t80, 0, -0.2e1 * t67, 0, 0, t49 * t68, t50 * t68, 0, 0, 0.2e1 * t76, -0.2e1 * t35 * t56 - 0.2e1 * t38 * t75, 0, 0.2e1 * t53, 0, 0, 0.2e1 * t40 * t75 + 0.2e1 * t56 * t46, 0.2e1 * t40 * t35 + 0.2e1 * t38 * t46, -0.2e1 * t12 * t38 - 0.2e1 * t13 * t56 - 0.2e1 * t15 * t35 - 0.2e1 * t16 * t75, 0.2e1 * t15 * t12 + 0.2e1 * t16 * t13 + 0.2e1 * t40 * t46, -0.2e1 * t78, 0.2e1 * t17 * t7 - 0.2e1 * t18 * t8, 0, 0.2e1 * t79, 0, 0, 0.2e1 * t22 * t17 + 0.2e1 * t21 * t8, 0.2e1 * t22 * t18 - 0.2e1 * t21 * t7, 0.2e1 * t1 * t17 - 0.2e1 * t2 * t18 + 0.2e1 * t3 * t7 - 0.2e1 * t4 * t8, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t21 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t56 + t13 * t38 - t15 * t75 + t16 * t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t18 - t2 * t17 - t3 * t8 - t4 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t53 + 0.2e1 * t76, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t78 + 0.2e1 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, 0, -t70, 0, -t44 * t69, t44 * t70, 0, 0, 0, 0, t35, 0, -t75, 0, t12, -t13, (-t73 * t35 - t72 * t75) * pkin(3), (t73 * t12 + t72 * t13) * pkin(3), 0, 0, -t7, 0, -t8, 0, t2, t1, -t23 * t17 + t24 * t18 + t32 * t7 - t33 * t8, -t1 * t33 + t2 * t32 + t4 * t23 - t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, -t69, 0, 0, 0, 0, 0, 0, 0, 0, -t75, -t35, 0, (t35 * t72 - t75 * t73) * pkin(3), 0, 0, 0, 0, 0, 0, -t8, t7, 0, t17 * t24 + t18 * t23 - t8 * t32 - t7 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t24, -0.2e1 * t23, 0, 0.2e1 * t33 * t23 - 0.2e1 * t32 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t35, 0, t46, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, -t23, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
