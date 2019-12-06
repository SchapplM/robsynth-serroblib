% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRPRR1
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
% Datum: 2019-12-05 15:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRPRR1_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR1_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRPRR1_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:43:04
% EndTime: 2019-12-05 15:43:07
% DurationCPUTime: 0.52s
% Computational Cost: add. (832->76), mult. (1930->148), div. (0->0), fcn. (1913->6), ass. (0->49)
t41 = cos(pkin(9));
t40 = sin(pkin(9));
t66 = cos(qJ(4));
t57 = t66 * t40;
t64 = sin(qJ(4));
t32 = t64 * t41 + t57;
t42 = sin(qJ(5));
t48 = t64 * t40 - t66 * t41;
t65 = cos(qJ(5));
t46 = t65 * t48;
t15 = t42 * t32 + t46;
t16 = t65 * t32 - t42 * t48;
t27 = t48 * qJD(4);
t54 = qJD(4) * t64;
t55 = qJD(4) * t66;
t61 = t40 * t55 + t41 * t54;
t8 = t16 * qJD(5) - t42 * t27 + t65 * t61;
t68 = t15 * t8;
t60 = qJD(5) * t42;
t7 = qJD(5) * t46 + t65 * t27 + t32 * t60 + t42 * t61;
t67 = t16 * t7;
t63 = t32 * t27;
t62 = pkin(6) + qJ(3);
t33 = t62 * t41;
t50 = t62 * t64;
t20 = t66 * t33 - t40 * t50;
t59 = pkin(4) * t60;
t37 = -t41 * pkin(3) - pkin(2);
t58 = t61 * pkin(4);
t56 = qJD(3) * t64;
t53 = t66 * qJD(3);
t52 = qJD(5) * t65 * pkin(4);
t51 = 0.2e1 * (t40 ^ 2 + t41 ^ 2) * qJD(3);
t30 = t62 * t57;
t19 = -t64 * t33 - t30;
t49 = t32 * pkin(7) - t19;
t47 = t65 * t49;
t45 = t48 * t61;
t12 = qJD(4) * t30 + t33 * t54 + t40 * t56 - t41 * t53;
t14 = -t48 * pkin(7) + t20;
t6 = t65 * t14 - t42 * t49;
t44 = t61 * pkin(7) + t12;
t13 = -t41 * t56 - t33 * t55 + (qJD(4) * t50 - t53) * t40;
t43 = t27 * pkin(7) + t13;
t21 = t48 * pkin(4) + t37;
t5 = -t42 * t14 - t47;
t2 = -t6 * qJD(5) + t42 * t44 + t65 * t43;
t1 = qJD(5) * t47 + t14 * t60 - t42 * t43 + t65 * t44;
t3 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t45 - 0.2e1 * t63, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t67 + 0.2e1 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32 * t12 - t48 * t13 - t61 * t19 - t27 * t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16 * t1 - t15 * t2 - t8 * t5 - t7 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, qJ(3) * t51, -0.2e1 * t63, 0.2e1 * t48 * t27 - 0.2e1 * t32 * t61, 0, 0.2e1 * t45, 0, 0, 0.2e1 * t37 * t61, -0.2e1 * t37 * t27, 0.2e1 * t12 * t48 - 0.2e1 * t13 * t32 + 0.2e1 * t19 * t27 - 0.2e1 * t20 * t61, -0.2e1 * t20 * t12 + 0.2e1 * t19 * t13, -0.2e1 * t67, 0.2e1 * t15 * t7 - 0.2e1 * t16 * t8, 0, 0.2e1 * t68, 0, 0, 0.2e1 * t15 * t58 + 0.2e1 * t21 * t8, 0.2e1 * t16 * t58 - 0.2e1 * t21 * t7, 0.2e1 * t1 * t15 - 0.2e1 * t2 * t16 + 0.2e1 * t5 * t7 - 0.2e1 * t6 * t8, -0.2e1 * t6 * t1 + 0.2e1 * t5 * t2 + 0.2e1 * t21 * t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, -t27, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t58; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, t27, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, (-t65 * t8 - t42 * t7 + (t15 * t42 + t65 * t16) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, 0, -t61, 0, t13, t12, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, (t65 * t7 - t42 * t8 + (-t65 * t15 + t16 * t42) * qJD(5)) * pkin(4), (t65 * t2 - t1 * t42 + (-t42 * t5 + t65 * t6) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t59, -0.2e1 * t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, -t8, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, -t52, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t3;
