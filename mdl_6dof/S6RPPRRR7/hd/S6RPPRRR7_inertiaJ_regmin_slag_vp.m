% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR7_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR7_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(pkin(10));
t40 = cos(pkin(10));
t44 = sin(qJ(4));
t46 = cos(qJ(4));
t24 = t39 * t46 + t40 * t44;
t25 = -t39 * t44 + t40 * t46;
t43 = sin(qJ(5));
t61 = cos(qJ(5));
t15 = t24 * t43 - t25 * t61;
t14 = t15 ^ 2;
t42 = sin(qJ(6));
t71 = t15 * t42;
t45 = cos(qJ(6));
t70 = t15 * t45;
t48 = t24 * t61 + t25 * t43;
t69 = t48 ^ 2;
t11 = t42 * t48;
t12 = t45 * t48;
t30 = pkin(3) * t39 + qJ(2);
t19 = pkin(4) * t24 + t30;
t68 = 0.2e1 * t19;
t67 = 0.2e1 * t30;
t66 = 0.2e1 * qJ(2);
t41 = -pkin(1) - qJ(3);
t62 = -pkin(7) + t41;
t26 = t62 * t39;
t27 = t62 * t40;
t54 = -t26 * t44 + t27 * t46;
t49 = -pkin(8) * t25 + t54;
t50 = -t26 * t46 - t27 * t44;
t9 = -pkin(8) * t24 - t50;
t4 = t43 * t9 - t49 * t61;
t65 = t4 * t45;
t64 = t43 * pkin(4);
t55 = t61 * pkin(4);
t32 = -t55 - pkin(5);
t63 = pkin(5) - t32;
t58 = t42 * t45;
t28 = t39 ^ 2 + t40 ^ 2;
t56 = 0.2e1 * t15 * t48;
t53 = pkin(5) * t15 - pkin(9) * t48;
t52 = -t69 - t14;
t31 = pkin(9) + t64;
t51 = -t15 * t32 - t31 * t48;
t47 = qJ(2) ^ 2;
t38 = t45 ^ 2;
t37 = t42 ^ 2;
t29 = 0.2e1 * t58;
t23 = t28 * t41;
t10 = t42 * t70;
t7 = (-t37 + t38) * t15;
t6 = pkin(5) * t48 + pkin(9) * t15 + t19;
t5 = t43 * t49 + t61 * t9;
t3 = t4 * t42;
t2 = t42 * t6 + t45 * t5;
t1 = -t42 * t5 + t45 * t6;
t8 = [1, 0, 0, -2 * pkin(1), t66 (pkin(1) ^ 2) + t47, t39 * t66, t40 * t66, -0.2e1 * t23, t28 * t41 ^ 2 + t47, t25 ^ 2, -0.2e1 * t25 * t24, 0, 0, 0, t24 * t67, t25 * t67, t14, t56, 0, 0, 0, t48 * t68, -t15 * t68, t38 * t14, -0.2e1 * t14 * t58, -0.2e1 * t48 * t70, t42 * t56, t69, 0.2e1 * t1 * t48 - 0.2e1 * t4 * t71, -0.2e1 * t2 * t48 - 0.2e1 * t4 * t70; 0, 0, 0, 1, 0, -pkin(1), 0, 0, -t28, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t52 * t42, t52 * t45; 0, 0, 0, 0, 0, 1, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t39, t40, 0, qJ(2), 0, 0, 0, 0, 0, t24, t25, 0, 0, 0, 0, 0, t48, -t15, 0, 0, 0, 0, 0, t12, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, t54, t50, 0, 0, -t15, -t48, 0, -t4, -t5, -t10, -t7, t11, t12, 0, t42 * t51 - t65, t45 * t51 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t24, 0, 0, 0, 0, 0, -t15, -t48, 0, 0, 0, 0, 0, -t70, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t55, -0.2e1 * t64, t37, t29, 0, 0, 0, -0.2e1 * t32 * t45, 0.2e1 * t32 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t48, 0, -t4, -t5, -t10, -t7, t11, t12, 0, t42 * t53 - t65, t45 * t53 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t48, 0, 0, 0, 0, 0, -t70, t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t55, -t64, t37, t29, 0, 0, 0, t63 * t45, -t63 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t37, t29, 0, 0, 0, 0.2e1 * pkin(5) * t45, -0.2e1 * pkin(5) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t70, t71, t48, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t42 * t31, -t45 * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t45, 0, -t42 * pkin(9), -t45 * pkin(9); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
