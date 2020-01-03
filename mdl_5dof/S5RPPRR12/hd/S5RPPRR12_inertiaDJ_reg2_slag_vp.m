% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR12
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
% Datum: 2019-12-31 18:07
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR12_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR12_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR12_inertiaDJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:07:21
% EndTime: 2019-12-31 18:07:24
% DurationCPUTime: 0.81s
% Computational Cost: add. (1025->100), mult. (2033->195), div. (0->0), fcn. (1898->6), ass. (0->73)
t31 = sin(pkin(8));
t33 = -pkin(1) - qJ(3);
t82 = -pkin(6) + t33;
t21 = t82 * t31;
t32 = cos(pkin(8));
t81 = sin(qJ(4));
t59 = t81 * qJD(3);
t60 = qJD(4) * t81;
t36 = cos(qJ(4));
t62 = t36 * t82;
t69 = t36 * qJD(3);
t88 = (qJD(4) * t62 - t59) * t32 - t21 * t60 - t31 * t69;
t19 = t36 * t31 + t81 * t32;
t20 = -t81 * t31 + t32 * t36;
t25 = pkin(3) * t31 + qJ(2);
t44 = pkin(4) * t19 - pkin(7) * t20 + t25;
t89 = -qJD(5) * t44 - t88;
t34 = sin(qJ(5));
t29 = t34 ^ 2;
t35 = cos(qJ(5));
t30 = t35 ^ 2;
t74 = t29 - t30;
t58 = qJD(5) * t74;
t22 = (t31 ^ 2 + t32 ^ 2) * qJD(3);
t56 = t82 * t81;
t13 = t36 * t21 + t32 * t56;
t15 = t19 * qJD(4);
t72 = qJD(4) * t36;
t16 = -t31 * t60 + t32 * t72;
t48 = pkin(4) * t16 + pkin(7) * t15 + qJD(2);
t71 = qJD(5) * t34;
t1 = t13 * t71 - t34 * t48 + t89 * t35;
t70 = qJD(5) * t35;
t2 = -t13 * t70 + t89 * t34 + t35 * t48;
t3 = -t34 * t13 + t35 * t44;
t4 = t35 * t13 + t34 * t44;
t51 = t3 * t34 - t35 * t4;
t87 = t51 * qJD(5) + t1 * t34 - t2 * t35;
t18 = t20 ^ 2;
t86 = 0.2e1 * qJD(2);
t85 = pkin(4) * t15;
t12 = t81 * t21 - t32 * t62;
t7 = t21 * t72 - t31 * t59 + (qJD(4) * t56 + t69) * t32;
t84 = t12 * t7;
t83 = t7 * t20;
t80 = t15 * t35;
t79 = t19 * t16;
t78 = t20 * t15;
t77 = t20 * t34;
t76 = t35 * t16;
t73 = t29 + t30;
t68 = qJ(2) * qJD(2);
t67 = 0.2e1 * t79;
t66 = -0.2e1 * pkin(4) * qJD(5);
t65 = t34 * t80;
t64 = t34 * t70;
t63 = t19 ^ 2 + t18;
t61 = t73 * t16;
t57 = t18 * t64;
t55 = -pkin(7) * t16 + t85;
t54 = pkin(4) * t20 + pkin(7) * t19;
t52 = t3 * t35 + t34 * t4;
t50 = t12 * t15 - t83;
t49 = t78 - t79;
t47 = -t15 * t34 + t20 * t70;
t46 = -t20 * t71 - t80;
t11 = t16 * t34 + t19 * t70;
t45 = 0.2e1 * t49;
t39 = -t52 * qJD(5) - t1 * t35 - t2 * t34;
t37 = t13 * t16 + t88 * t19 + t50;
t10 = t19 * t71 - t76;
t5 = t20 * t58 + t65;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, 0.2e1 * t68, 0, 0, 0, 0, 0, 0, t31 * t86, t32 * t86, 0.2e1 * t22, -0.2e1 * t33 * t22 + 0.2e1 * t68, -0.2e1 * t78, 0.2e1 * t15 * t19 - 0.2e1 * t16 * t20, 0, t67, 0, 0, 0.2e1 * qJD(2) * t19 + 0.2e1 * t16 * t25, 0.2e1 * qJD(2) * t20 - 0.2e1 * t15 * t25, -0.2e1 * t37, 0.2e1 * t25 * qJD(2) + 0.2e1 * t88 * t13 + 0.2e1 * t84, -0.2e1 * t30 * t78 - 0.2e1 * t57, 0.2e1 * t18 * t58 + 0.4e1 * t20 * t65, 0.2e1 * t46 * t19 + 0.2e1 * t20 * t76, -0.2e1 * t29 * t78 + 0.2e1 * t57, -0.2e1 * t16 * t77 - 0.2e1 * t47 * t19, t67, 0.2e1 * t47 * t12 + 0.2e1 * t16 * t3 + 0.2e1 * t19 * t2 + 0.2e1 * t7 * t77, 0.2e1 * t1 * t19 + 0.2e1 * t46 * t12 - 0.2e1 * t16 * t4 + 0.2e1 * t35 * t83, 0.2e1 * t52 * t15 + 0.2e1 * t87 * t20, -0.2e1 * t1 * t4 + 0.2e1 * t2 * t3 + 0.2e1 * t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, t45, t37, 0, 0, 0, 0, 0, 0, t34 * t45 - t63 * t70, t35 * t45 + t63 * t71, 0, -t51 * t16 + t39 * t19 + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t19 * t61 - 0.2e1 * t78; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(2), 0, 0, 0, 0, 0, 0, t16, -t15, 0, qJD(2), 0, 0, 0, 0, 0, 0, -t10, -t11, t73 * t15, -t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, 0, -t16, 0, -t7, -t88, 0, 0, -t5, t74 * t15 - 0.4e1 * t20 * t64, t11, t5, -t10, 0, -t35 * t7 + t55 * t34 + (t12 * t34 - t54 * t35) * qJD(5), t34 * t7 + t55 * t35 + (t12 * t35 + t54 * t34) * qJD(5), t39, -pkin(4) * t7 + t39 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t47, t61, pkin(7) * t61 - t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t64, -0.2e1 * t58, 0, -0.2e1 * t64, 0, 0, t34 * t66, t35 * t66, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, 0, -t47, t16, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, -t70, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, 0, -t71, 0, -pkin(7) * t70, pkin(7) * t71, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t6;
