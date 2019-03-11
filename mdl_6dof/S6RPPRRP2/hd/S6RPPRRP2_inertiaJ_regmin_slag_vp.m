% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRP2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = sin(pkin(10));
t36 = cos(pkin(10));
t39 = sin(qJ(4));
t58 = cos(qJ(4));
t23 = t34 * t58 + t36 * t39;
t67 = -0.2e1 * t23;
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t43 = pkin(5) * t38 - qJ(6) * t40;
t66 = t43 * t23;
t37 = cos(pkin(9));
t29 = -pkin(1) * t37 - pkin(2);
t24 = -pkin(3) * t36 + t29;
t65 = 0.2e1 * t24;
t64 = -0.2e1 * t38;
t35 = sin(pkin(9));
t27 = pkin(1) * t35 + qJ(3);
t59 = pkin(7) + t27;
t18 = t59 * t34;
t19 = t59 * t36;
t11 = -t18 * t39 + t19 * t58;
t22 = t34 * t39 - t36 * t58;
t9 = pkin(4) * t22 - pkin(8) * t23 + t24;
t4 = t11 * t40 + t38 * t9;
t63 = pkin(8) * t22;
t62 = t22 * pkin(5);
t61 = t38 * pkin(8);
t60 = t40 * pkin(8);
t32 = t38 ^ 2;
t57 = t32 * t23;
t13 = t38 * t22;
t56 = t38 * t23;
t55 = t38 * t40;
t16 = t40 * t22;
t17 = t40 * t23;
t54 = t34 ^ 2 + t36 ^ 2;
t33 = t40 ^ 2;
t53 = t32 + t33;
t52 = t22 * qJ(6);
t51 = t22 * t67;
t50 = t53 * pkin(8);
t49 = t11 * t38 - t40 * t9;
t48 = -pkin(4) * t23 - t63;
t1 = t52 + t4;
t2 = t49 - t62;
t47 = t1 * t40 + t2 * t38;
t46 = t1 * t38 - t2 * t40;
t44 = pkin(5) * t40 + qJ(6) * t38;
t25 = -pkin(4) - t44;
t45 = -t23 * t25 + t63;
t10 = t18 * t58 + t39 * t19;
t21 = t23 ^ 2;
t20 = t22 ^ 2;
t15 = t33 * t23;
t14 = t33 * t21;
t12 = -t15 - t57;
t5 = t10 + t66;
t3 = [1, 0, 0 (t35 ^ 2 + t37 ^ 2) * pkin(1) ^ 2, -0.2e1 * t29 * t36, 0.2e1 * t29 * t34, 0.2e1 * t54 * t27, t27 ^ 2 * t54 + t29 ^ 2, t21, t51, 0, 0, 0, t22 * t65, t23 * t65, t14, -0.2e1 * t21 * t55, 0.2e1 * t22 * t17, t38 * t51, t20, 0.2e1 * t10 * t56 - 0.2e1 * t22 * t49, 0.2e1 * t10 * t17 - 0.2e1 * t22 * t4, -0.2e1 * t2 * t22 + 0.2e1 * t5 * t56, t46 * t67, 0.2e1 * t1 * t22 - 0.2e1 * t17 * t5, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t22 + t23 * t47; 0, 0, 0, 1, 0, 0, 0, t54, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t32 + t14 + t20; 0, 0, 0, 0, -t36, t34, 0, t29, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, t16, -t13, t16, t12, t13, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, -t10, -t11, t38 * t17, t15 - t57, t13, t16, 0, -t10 * t40 + t38 * t48, t10 * t38 + t40 * t48, -t38 * t45 - t5 * t40, t47, -t5 * t38 + t40 * t45, pkin(8) * t47 + t5 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, 0, 0, 0, 0, 0, -t16, t13, -t16, -t12, -t13, t22 * t25 + t23 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t32, 0.2e1 * t55, 0, 0, 0, 0.2e1 * pkin(4) * t40, pkin(4) * t64, -0.2e1 * t25 * t40, 0.2e1 * t50, t25 * t64, pkin(8) ^ 2 * t53 + t25 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t56, t22, -t49, -t4, -t49 + 0.2e1 * t62, -t44 * t23, 0.2e1 * t52 + t4, -pkin(5) * t2 + qJ(6) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t56, -t17, -t56, 0, t17, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, t40, 0, t38, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t40, 0, -t61, -t60, -t61, -t43, t60, -t43 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, t17, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t3;
