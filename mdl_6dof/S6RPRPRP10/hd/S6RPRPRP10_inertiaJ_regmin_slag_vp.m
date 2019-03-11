% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:32
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP10_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPRP10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(3));
t41 = cos(qJ(3));
t17 = t39 * pkin(3) - t41 * qJ(4) + qJ(2);
t58 = -0.2e1 * t17;
t57 = 0.2e1 * qJ(2);
t56 = 0.2e1 * qJ(4);
t43 = -pkin(1) - pkin(7);
t19 = (pkin(4) - t43) * t41;
t38 = sin(qJ(5));
t40 = cos(qJ(5));
t8 = t39 * pkin(8) + t17;
t5 = t38 * t19 + t40 * t8;
t55 = t41 * pkin(5);
t25 = t38 * t39;
t42 = -pkin(3) - pkin(8);
t54 = t38 * t42;
t44 = -t38 * pkin(5) + t40 * qJ(6);
t16 = qJ(4) - t44;
t53 = t39 * t16;
t52 = t40 * t38;
t27 = t40 * t39;
t51 = t41 * t39;
t50 = t41 * t42;
t49 = t41 * t43;
t32 = t38 ^ 2;
t34 = t40 ^ 2;
t23 = t32 + t34;
t33 = t39 ^ 2;
t35 = t41 ^ 2;
t24 = t33 + t35;
t48 = t39 * qJ(4);
t47 = t41 * qJ(6);
t46 = 0.2e1 * t51;
t45 = -t40 * t19 + t38 * t8;
t2 = t47 + t5;
t3 = t45 - t55;
t1 = t2 * t38 - t3 * t40;
t21 = t41 * pkin(3) + t48;
t20 = pkin(5) * t40 + t38 * qJ(6);
t30 = t40 * t42;
t29 = t39 * t43;
t28 = t40 * t41;
t26 = t38 * t41;
t22 = t40 * t50;
t18 = -t39 * pkin(4) + t29;
t15 = t24 * t43;
t14 = t23 * t42;
t13 = t24 * t40;
t12 = t24 * t38;
t11 = t23 * t41;
t6 = t29 + (-pkin(4) - t20) * t39;
t4 = [1, 0, 0, -2 * pkin(1), t57 (pkin(1) ^ 2) + qJ(2) ^ 2, t35, -0.2e1 * t51, 0, 0, 0, t39 * t57, t41 * t57, -0.2e1 * t15, t39 * t58, t41 * t58, t24 * t43 ^ 2 + t17 ^ 2, t32 * t33, 0.2e1 * t33 * t52, t38 * t46, t40 * t46, t35, -0.2e1 * t18 * t27 - 0.2e1 * t41 * t45, 0.2e1 * t18 * t25 - 0.2e1 * t5 * t41, -0.2e1 * t6 * t27 - 0.2e1 * t3 * t41, 0.2e1 * (t2 * t40 + t3 * t38) * t39, 0.2e1 * t2 * t41 - 0.2e1 * t6 * t25, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, -t24, 0, 0, t15, 0, 0, 0, 0, 0, -t13, t12, -t13, 0, -t12, -t1 * t41 + t6 * t39; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t35 + t33; 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, t49, -t29, -t21, -t49, t29, t21 * t43, t38 * t27 (-t32 + t34) * t39, t28, -t26, 0, t18 * t38 - t40 * t48 + t22, t18 * t40 + (t48 - t50) * t38, -t16 * t27 + t6 * t38 + t22, -t1, -t6 * t40 + (t50 - t53) * t38, t1 * t42 + t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, -t39, 0, -t41, t39, t21, 0, 0, 0, 0, 0, t25, t27, t25, t11, -t27, -t23 * t50 + t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t56, pkin(3) ^ 2 + qJ(4) ^ 2, t34, -0.2e1 * t52, 0, 0, 0, t38 * t56, t40 * t56, 0.2e1 * t16 * t38, -0.2e1 * t14, -0.2e1 * t16 * t40, t23 * t42 ^ 2 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, -t49, 0, 0, 0, 0, 0, t28, -t26, t28, 0, t26, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, -t23, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t27, t41, -t45, -t5, -t45 + 0.2e1 * t55, t44 * t39, 0.2e1 * t47 + t5, -t3 * pkin(5) + t2 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t26, -t28, 0, -t26, -t20 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, t30, -t54, t30, -t20, t54, t20 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, t40, 0, t38, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(5), 0, 0.2e1 * qJ(6), pkin(5) ^ 2 + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, t25, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, 0, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t4;
