% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S5RRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RRRPR2_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR2_inertiaDJ_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR2_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:07:37
% EndTime: 2020-01-03 12:07:38
% DurationCPUTime: 0.25s
% Computational Cost: add. (331->62), mult. (884->102), div. (0->0), fcn. (618->8), ass. (0->54)
t46 = cos(qJ(2));
t37 = t46 * pkin(1) + pkin(2);
t45 = cos(qJ(3));
t59 = pkin(1) * qJD(2);
t52 = t46 * t59;
t42 = sin(qJ(3));
t43 = sin(qJ(2));
t54 = t42 * t43 * pkin(1);
t56 = qJD(3) * t45;
t11 = -t37 * t56 - t45 * t52 + (qJD(2) + qJD(3)) * t54;
t61 = t43 * t45;
t47 = (-t43 * t56 + (-t42 * t46 - t61) * qJD(2)) * pkin(1);
t57 = qJD(3) * t42;
t12 = -t37 * t57 + t47;
t39 = sin(pkin(9));
t40 = cos(pkin(9));
t2 = -t39 * t11 - t40 * t12;
t44 = cos(qJ(5));
t38 = t44 * qJD(5);
t41 = sin(qJ(5));
t22 = t45 * t37 + pkin(3) - t54;
t25 = pkin(1) * t61 + t42 * t37;
t8 = t40 * t22 - t39 * t25;
t6 = -pkin(4) - t8;
t64 = t2 * t41 + t6 * t38;
t63 = t39 * t42;
t62 = t40 * t42;
t36 = t45 * pkin(2) + pkin(3);
t23 = -pkin(2) * t63 + t40 * t36;
t18 = -pkin(4) - t23;
t58 = pkin(2) * qJD(3);
t20 = (t39 * t45 + t62) * t58;
t60 = t18 * t38 + t20 * t41;
t9 = t39 * t22 + t40 * t25;
t24 = pkin(2) * t62 + t39 * t36;
t55 = t41 * qJD(5);
t53 = t43 * t59;
t51 = pkin(2) * t57;
t50 = pkin(2) * t56;
t4 = t6 * t55;
t49 = -t2 * t44 + t4;
t13 = t18 * t55;
t48 = -t20 * t44 + t13;
t35 = -t40 * pkin(3) - pkin(4);
t34 = t39 * pkin(3) + pkin(8);
t32 = 0.2e1 * t41 * t38;
t28 = t35 * t38;
t27 = t35 * t55;
t26 = 0.2e1 * (-t41 ^ 2 + t44 ^ 2) * qJD(5);
t21 = (t40 * t45 - t63) * t58;
t19 = pkin(8) + t24;
t7 = pkin(8) + t9;
t3 = -t40 * t11 + t39 * t12;
t1 = [0, 0, 0, 0, -0.2e1 * t53, -0.2e1 * t52, 0, 0.2e1 * t12, 0.2e1 * t11, -0.2e1 * t8 * t2 + 0.2e1 * t9 * t3, t32, t26, 0, 0, 0, 0.2e1 * t49, 0.2e1 * t64; 0, 0, 0, 0, -t53, -t52, 0, (-pkin(2) - t37) * t57 + t47, t11 - t50, -t2 * t23 - t8 * t20 + t9 * t21 + t3 * t24, t32, t26, 0, 0, 0, t13 + t4 + (-t2 - t20) * t44, t60 + t64; 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t51, -0.2e1 * t50, -0.2e1 * t23 * t20 + 0.2e1 * t24 * t21, t32, t26, 0, 0, 0, 0.2e1 * t48, 0.2e1 * t60; 0, 0, 0, 0, 0, 0, 0, t12, t11, (-t2 * t40 + t3 * t39) * pkin(3), t32, t26, 0, 0, 0, t27 + t49, t28 + t64; 0, 0, 0, 0, 0, 0, 0, -t51, -t50, (-t20 * t40 + t21 * t39) * pkin(3), t32, t26, 0, 0, 0, t27 + t48, t28 + t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t26, 0, 0, 0, 0.2e1 * t27, 0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t55, 0, -t41 * t3 - t7 * t38, -t44 * t3 + t7 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t55, 0, -t19 * t38 - t41 * t21, t19 * t55 - t44 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t55, 0, -t34 * t38, t34 * t55; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t1;
