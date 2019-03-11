% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(pkin(9));
t38 = sin(pkin(10));
t40 = cos(pkin(10));
t49 = t38 ^ 2 + t40 ^ 2;
t56 = t49 * t39;
t43 = sin(qJ(5));
t53 = cos(qJ(5));
t19 = t53 * t38 + t43 * t40;
t55 = -0.2e1 * t19;
t41 = cos(pkin(9));
t45 = -pkin(1) - pkin(2);
t25 = t41 * qJ(2) + t39 * t45;
t21 = -qJ(4) + t25;
t54 = pkin(7) - t21;
t44 = cos(qJ(6));
t46 = -t43 * t38 + t53 * t40;
t14 = t46 * t44;
t42 = sin(qJ(6));
t13 = t42 * t46;
t52 = t42 * t19;
t51 = t42 * t44;
t50 = t44 * t19;
t23 = t39 * qJ(2) - t41 * t45;
t48 = t46 * t55;
t22 = pkin(3) + t23;
t12 = t40 * pkin(4) + t22;
t47 = pkin(5) * t19 - pkin(8) * t46;
t37 = t44 ^ 2;
t36 = t42 ^ 2;
t35 = t41 ^ 2;
t33 = t39 ^ 2;
t15 = t19 ^ 2;
t11 = t46 * t39;
t10 = t19 * t39;
t9 = t54 * t40;
t8 = t54 * t38;
t7 = t44 * t11 - t42 * t41;
t6 = -t42 * t11 - t44 * t41;
t5 = pkin(5) * t46 + t19 * pkin(8) + t12;
t4 = t43 * t8 - t53 * t9;
t3 = -t43 * t9 - t53 * t8;
t2 = t44 * t4 + t42 * t5;
t1 = -t42 * t4 + t44 * t5;
t16 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2) (pkin(1) ^ 2) + qJ(2) ^ 2, 0.2e1 * t23, 0.2e1 * t25, t23 ^ 2 + t25 ^ 2, 0.2e1 * t22 * t40, -0.2e1 * t22 * t38, -0.2e1 * t49 * t21, t49 * t21 ^ 2 + t22 ^ 2, t15, -t48, 0, 0, 0, 0.2e1 * t12 * t46, t12 * t55, t37 * t15, -0.2e1 * t15 * t51, -0.2e1 * t46 * t50, -t42 * t48, t46 ^ 2, 0.2e1 * t1 * t46 - 0.2e1 * t3 * t52, -0.2e1 * t2 * t46 - 0.2e1 * t3 * t50; 0, 0, 0, -1, 0, -pkin(1), -t41, t39, -t23 * t41 + t25 * t39, -t41 * t40, t41 * t38, -t56, t21 * t56 - t22 * t41, 0, 0, 0, 0, 0, -t41 * t46, t41 * t19, 0, 0, 0, 0, 0, -t10 * t52 + t46 * t6, -t10 * t50 - t46 * t7; 0, 0, 0, 0, 0, 1, 0, 0, t33 + t35, 0, 0, 0, t49 * t33 + t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t38, 0, t22, 0, 0, 0, 0, 0, t46, -t19, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19, -t46, 0, -t3, -t4, -t42 * t50 (t36 - t37) * t19, t13, t14, 0, -t3 * t44 + t47 * t42, t3 * t42 + t47 * t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t10 * t44, t10 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t19, 0, 0, 0, 0, 0, t14, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t36, 0.2e1 * t51, 0, 0, 0, 0.2e1 * pkin(5) * t44, -0.2e1 * pkin(5) * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, t52, t46, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t52, -t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t44, 0, -t42 * pkin(8), -t44 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t16;
