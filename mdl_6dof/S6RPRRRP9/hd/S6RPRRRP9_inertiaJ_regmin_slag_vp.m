% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x29]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:29
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRRRP9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP9_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRRP9_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:53:24
% EndTime: 2019-05-06 01:53:25
% DurationCPUTime: 0.64s
% Computational Cost: add. (546->104), mult. (1069->186), div. (0->0), fcn. (1174->6), ass. (0->61)
t40 = sin(qJ(5));
t41 = sin(qJ(4));
t43 = cos(qJ(5));
t44 = cos(qJ(4));
t23 = t40 * t44 + t43 * t41;
t45 = cos(qJ(3));
t17 = t45 * t23;
t67 = -0.2e1 * t17;
t34 = -t44 * pkin(4) - pkin(3);
t66 = 0.2e1 * t34;
t65 = 0.2e1 * t44;
t64 = 2 * qJ(2);
t63 = pkin(8) + pkin(9);
t62 = t40 * pkin(4);
t35 = t43 * pkin(4);
t42 = sin(qJ(3));
t61 = t23 * t42;
t60 = t41 * t42;
t59 = t41 * t44;
t58 = t41 * t45;
t46 = -pkin(1) - pkin(7);
t57 = t41 * t46;
t56 = t42 * t46;
t25 = t42 * pkin(3) - t45 * pkin(8) + qJ(2);
t53 = t44 * t46;
t48 = t42 * t53;
t10 = t48 + (-pkin(9) * t45 + t25) * t41;
t55 = t43 * t10;
t54 = t44 * t42;
t32 = t44 * t45;
t52 = t45 * t42;
t51 = t45 * t46;
t37 = t42 ^ 2;
t39 = t45 ^ 2;
t50 = -t37 - t39;
t49 = -0.2e1 * t52;
t20 = t44 * t25;
t8 = -pkin(9) * t32 + t20 + (pkin(4) - t57) * t42;
t3 = -t40 * t10 + t43 * t8;
t26 = t63 * t41;
t27 = t63 * t44;
t11 = -t43 * t26 - t40 * t27;
t21 = pkin(4) * t58 - t51;
t47 = -pkin(3) * t45 - pkin(8) * t42;
t4 = t40 * t8 + t55;
t12 = -t40 * t26 + t43 * t27;
t38 = t44 ^ 2;
t36 = t41 ^ 2;
t33 = t35 + pkin(5);
t22 = t40 * t41 - t43 * t44;
t19 = t43 * t32 - t40 * t58;
t18 = -t40 * t60 + t43 * t54;
t15 = t22 * pkin(5) + t34;
t14 = t41 * t25 + t48;
t13 = -t41 * t56 + t20;
t9 = t17 * pkin(5) + t21;
t6 = -t22 * qJ(6) + t12;
t5 = -t23 * qJ(6) + t11;
t2 = -t17 * qJ(6) + t4;
t1 = t42 * pkin(5) - t19 * qJ(6) + t3;
t7 = [1, 0, 0, -2 * pkin(1), t64, pkin(1) ^ 2 + qJ(2) ^ 2, t39, t49, 0, 0, 0, t42 * t64, t45 * t64, t38 * t39, -0.2e1 * t39 * t59, t52 * t65, t41 * t49, t37, 0.2e1 * t13 * t42 - 0.2e1 * t39 * t57, -0.2e1 * t14 * t42 - 0.2e1 * t39 * t53, t19 ^ 2, t19 * t67, 0.2e1 * t19 * t42, t42 * t67, t37, 0.2e1 * t21 * t17 + 0.2e1 * t3 * t42, 0.2e1 * t21 * t19 - 0.2e1 * t4 * t42, -0.2e1 * t1 * t19 - 0.2e1 * t2 * t17, t1 ^ 2 + t2 ^ 2 + t9 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t50 * t41, t50 * t44, 0, 0, 0, 0, 0, -t45 * t17 - t42 * t61, -t18 * t42 - t45 * t19, -t18 * t17 + t19 * t61, -t1 * t61 + t2 * t18 - t9 * t45; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 ^ 2 + t61 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, t51, -t56, t41 * t32 (-t36 + t38) * t45, t60, t54, 0, t47 * t41 + t44 * t51, -t41 * t51 + t47 * t44, t19 * t23, -t23 * t17 - t19 * t22, t61, -t22 * t42, 0, t11 * t42 + t34 * t17 + t21 * t22, -t12 * t42 + t34 * t19 + t21 * t23, -t1 * t23 - t6 * t17 - t5 * t19 - t2 * t22, t1 * t5 + t9 * t15 + t2 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42, 0, 0, 0, 0, 0, t32, -t58, 0, 0, 0, 0, 0, -t45 * t22, -t17, -t18 * t22 + t23 * t61, -t45 * t15 + t18 * t6 - t5 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, 0.2e1 * t59, 0, 0, 0, pkin(3) * t65, -0.2e1 * pkin(3) * t41, t23 ^ 2, -0.2e1 * t23 * t22, 0, 0, 0, t22 * t66, t23 * t66, -0.2e1 * t6 * t22 - 0.2e1 * t5 * t23, t15 ^ 2 + t5 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t58, t42, t13, -t14, 0, 0, t19, -t17, t42, t42 * t35 + t3, -t55 + (-t42 * pkin(4) - t8) * t40, -t17 * t62 - t33 * t19, t1 * t33 + t2 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, -t54, 0, 0, 0, 0, 0, -t61, -t18, 0, t18 * t62 - t33 * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, t44, 0, -t41 * pkin(8), -t44 * pkin(8), 0, 0, t23, -t22, 0, t11, -t12, -t22 * t62 - t33 * t23, t5 * t33 + t6 * t62; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t62, 0, t40 ^ 2 * pkin(4) ^ 2 + t33 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, t42, t3, -t4, -pkin(5) * t19, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t61, -t18, 0, -t61 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t22, 0, t11, -t12, -pkin(5) * t23, t5 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t35, -t62, 0, t33 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t7;
