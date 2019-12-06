% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5PRRRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:06:18
% EndTime: 2019-12-05 17:06:21
% DurationCPUTime: 0.52s
% Computational Cost: add. (260->56), mult. (804->95), div. (0->0), fcn. (495->6), ass. (0->49)
t35 = cos(qJ(4));
t58 = cos(qJ(3));
t46 = t58 * pkin(2);
t39 = t46 + pkin(3);
t37 = t35 * t39;
t40 = qJD(3) * t46;
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t59 = pkin(2) * t33;
t52 = t32 * t59;
t61 = qJD(3) + qJD(4);
t5 = -qJD(4) * t37 - t35 * t40 + t52 * t61;
t31 = sin(qJ(5));
t29 = t31 ^ 2;
t34 = cos(qJ(5));
t30 = t34 ^ 2;
t55 = t29 + t30;
t43 = t55 * t5;
t57 = t33 * t35;
t63 = t61 * pkin(2) * (t32 * t58 + t57);
t28 = t34 * qJD(5);
t54 = qJD(4) * pkin(3);
t48 = t32 * t54;
t6 = t48 + t63;
t12 = t37 - t52;
t9 = -pkin(4) - t12;
t60 = t9 * t28 + t6 * t31;
t47 = t35 * t54;
t11 = t55 * t47;
t27 = -t35 * pkin(3) - pkin(4);
t56 = t27 * t28 + t31 * t48;
t13 = pkin(2) * t57 + t32 * t39;
t53 = t31 * qJD(5);
t51 = pkin(4) * t53;
t50 = pkin(4) * t28;
t49 = qJD(3) * t59;
t45 = t31 * t28;
t7 = t9 * t53;
t44 = -t34 * t6 + t7;
t42 = t55 * t35;
t41 = -0.2e1 * t48;
t15 = t27 * t53;
t38 = -t34 * t48 + t15;
t26 = t32 * pkin(3) + pkin(8);
t24 = -0.2e1 * t45;
t23 = 0.2e1 * t45;
t14 = 0.2e1 * (-t29 + t30) * qJD(5);
t10 = pkin(8) + t13;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t49, -0.2e1 * t40, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t6, 0.2e1 * t5, 0, -0.2e1 * t12 * t6 - 0.2e1 * t13 * t5, t23, t14, 0, t24, 0, 0, 0.2e1 * t44, 0.2e1 * t60, -0.2e1 * t43, -0.2e1 * t10 * t43 + 0.2e1 * t9 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t41 - t63, t5 - t47, 0, (-t32 * t5 - t35 * t6 + (-t12 * t32 + t13 * t35) * qJD(4)) * pkin(3), t23, t14, 0, t24, 0, 0, t15 + t7 + (-t6 - t48) * t34, t56 + t60, t11 - t43, t6 * t27 - t26 * t43 + (t10 * t42 + t32 * t9) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -0.2e1 * t47, 0, 0, t23, t14, 0, t24, 0, 0, 0.2e1 * t38, 0.2e1 * t56, 0.2e1 * t11, 0.2e1 * (t26 * t42 + t27 * t32) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, 0, 0, t23, t14, 0, t24, 0, 0, t44 - t51, -t50 + t60, -t43, -t6 * pkin(4) - pkin(8) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, 0, 0, t23, t14, 0, t24, 0, 0, t38 - t51, -t50 + t56, t11, (-pkin(4) * t32 + pkin(8) * t42) * t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, t14, 0, t24, 0, 0, -0.2e1 * t51, -0.2e1 * t50, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t53, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -t10 * t28 + t31 * t5, t10 * t53 + t34 * t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -t26 * t28 - t31 * t47, t26 * t53 - t34 * t47, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, 0, -t53, 0, -pkin(8) * t28, pkin(8) * t53, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
