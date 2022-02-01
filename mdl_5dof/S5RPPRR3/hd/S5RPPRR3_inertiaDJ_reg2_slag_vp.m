% Calculate inertial parameters regressor of joint inertia matrix time derivative for
% S5RPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MMD_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of inertia matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:15
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S5RPPRR3_inertiaDJ_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR3_inertiaDJ_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR3_inertiaDJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:14:36
% EndTime: 2022-01-23 09:14:38
% DurationCPUTime: 0.50s
% Computational Cost: add. (995->78), mult. (2093->150), div. (0->0), fcn. (2076->8), ass. (0->50)
t42 = cos(pkin(9));
t41 = sin(pkin(9));
t68 = cos(qJ(4));
t59 = t68 * t41;
t66 = sin(qJ(4));
t33 = t66 * t42 + t59;
t44 = sin(qJ(5));
t50 = t66 * t41 - t68 * t42;
t67 = cos(qJ(5));
t18 = t67 * t33 - t44 * t50;
t31 = t50 * qJD(4);
t48 = t67 * t50;
t62 = qJD(5) * t44;
t56 = qJD(4) * t66;
t57 = qJD(4) * t68;
t63 = t41 * t57 + t42 * t56;
t9 = qJD(5) * t48 + t67 * t31 + t33 * t62 + t44 * t63;
t70 = t18 * t9;
t38 = sin(pkin(8)) * pkin(1) + qJ(3);
t69 = pkin(6) + t38;
t10 = t18 * qJD(5) - t44 * t31 + t67 * t63;
t17 = t44 * t33 + t48;
t65 = t17 * t10;
t64 = t33 * t31;
t32 = t69 * t42;
t52 = t69 * t66;
t16 = t68 * t32 - t41 * t52;
t61 = pkin(4) * t62;
t60 = t63 * pkin(4);
t58 = qJD(3) * t66;
t55 = t68 * qJD(3);
t54 = qJD(5) * t67 * pkin(4);
t53 = 0.2e1 * (t41 ^ 2 + t42 ^ 2) * qJD(3);
t27 = t69 * t59;
t15 = -t66 * t32 - t27;
t34 = -cos(pkin(8)) * pkin(1) - t42 * pkin(3) - pkin(2);
t51 = t33 * pkin(7) - t15;
t49 = t67 * t51;
t47 = t50 * t63;
t12 = qJD(4) * t27 + t32 * t56 + t41 * t58 - t42 * t55;
t14 = -t50 * pkin(7) + t16;
t4 = t67 * t14 - t44 * t51;
t46 = t63 * pkin(7) + t12;
t13 = -t42 * t58 - t32 * t57 + (qJD(4) * t52 - t55) * t41;
t45 = t31 * pkin(7) + t13;
t21 = t50 * pkin(4) + t34;
t3 = -t44 * t14 - t49;
t2 = -t4 * qJD(5) + t44 * t46 + t67 * t45;
t1 = qJD(5) * t49 + t14 * t62 - t44 * t45 + t67 * t46;
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t38 * t53, -0.2e1 * t64, 0.2e1 * t50 * t31 - 0.2e1 * t33 * t63, 0, 0.2e1 * t47, 0, 0, 0.2e1 * t34 * t63, -0.2e1 * t34 * t31, 0.2e1 * t12 * t50 - 0.2e1 * t13 * t33 + 0.2e1 * t15 * t31 - 0.2e1 * t16 * t63, -0.2e1 * t16 * t12 + 0.2e1 * t15 * t13, -0.2e1 * t70, -0.2e1 * t18 * t10 + 0.2e1 * t17 * t9, 0, 0.2e1 * t65, 0, 0, 0.2e1 * t21 * t10 + 0.2e1 * t17 * t60, 0.2e1 * t18 * t60 - 0.2e1 * t21 * t9, 0.2e1 * t1 * t17 - 0.2e1 * t4 * t10 - 0.2e1 * t2 * t18 + 0.2e1 * t3 * t9, -0.2e1 * t4 * t1 + 0.2e1 * t3 * t2 + 0.2e1 * t21 * t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12 * t33 - t13 * t50 - t15 * t63 - t16 * t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t18 - t3 * t10 - t2 * t17 - t4 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t47 - 0.2e1 * t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t65 - 0.2e1 * t70; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, -t31, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t60; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, 0, -t63, 0, t13, t12, 0, 0, 0, 0, -t9, 0, -t10, 0, t2, t1, (t67 * t9 - t10 * t44 + (-t67 * t17 + t18 * t44) * qJD(5)) * pkin(4), (t67 * t2 - t1 * t44 + (-t3 * t44 + t67 * t4) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t63, t31, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, (-t67 * t10 - t44 * t9 + (t17 * t44 + t67 * t18) * qJD(5)) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t61, -0.2e1 * t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, -t10, 0, t2, t1, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t54, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg = t5;
