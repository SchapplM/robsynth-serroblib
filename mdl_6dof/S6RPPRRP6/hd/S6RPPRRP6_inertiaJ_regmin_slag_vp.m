% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:11
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t57 = 2 * pkin(5);
t22 = (pkin(1) + qJ(3));
t56 = 2 * t22;
t24 = sin(qJ(5));
t55 = -0.2e1 * t24;
t26 = cos(qJ(5));
t54 = 0.2e1 * t26;
t27 = cos(qJ(4));
t53 = 0.2e1 * t27;
t25 = sin(qJ(4));
t52 = pkin(8) * t25;
t51 = t24 * pkin(8);
t50 = t26 * pkin(8);
t10 = t25 * pkin(4) - t27 * pkin(8) + t22;
t15 = t26 * t25;
t21 = -pkin(7) + qJ(2);
t4 = t24 * t10 + t21 * t15;
t17 = t24 ^ 2;
t49 = t17 * t27;
t48 = t21 * t24;
t47 = t24 * t26;
t46 = t24 * t27;
t45 = t25 * t21;
t16 = t26 * t27;
t32 = -t26 * pkin(5) - t24 * qJ(6);
t11 = -pkin(4) + t32;
t44 = t27 * t11;
t43 = t27 * t21;
t42 = t27 * t25;
t19 = t26 ^ 2;
t41 = t17 + t19;
t18 = t25 ^ 2;
t20 = t27 ^ 2;
t40 = t18 + t20;
t39 = t25 * qJ(6);
t38 = -0.2e1 * t42;
t37 = t41 * t25;
t36 = -pkin(4) * t27 - t52;
t1 = t39 + t4;
t7 = t26 * t10;
t2 = -t7 + (-pkin(5) + t48) * t25;
t35 = t1 * t26 + t2 * t24;
t34 = -t1 * t24 + t2 * t26;
t33 = -t44 + t52;
t31 = -pkin(5) * t24 + t26 * qJ(6);
t29 = (qJ(2) ^ 2);
t28 = 2 * qJ(2);
t14 = t19 * t27;
t13 = t24 * t25;
t9 = t40 * t26;
t8 = t40 * t24;
t5 = (-t21 - t31) * t27;
t3 = -t24 * t45 + t7;
t6 = [1, 0, 0, -2 * pkin(1), t28, pkin(1) ^ 2 + t29, t28, t56, t22 ^ 2 + t29, t20, t38, 0, 0, 0, t25 * t56, t22 * t53, t19 * t20, -0.2e1 * t20 * t47, t42 * t54, t24 * t38, t18, -0.2e1 * t20 * t48 + 0.2e1 * t3 * t25, -0.2e1 * t20 * t21 * t26 - 0.2e1 * t4 * t25, -0.2e1 * t2 * t25 + 0.2e1 * t5 * t46, t34 * t53, 0.2e1 * t1 * t25 - 0.2e1 * t5 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t22, 0, 0, 0, 0, 0, -t25, -t27, 0, 0, 0, 0, 0, -t15, t13, -t15, t14 + t49, -t13, t34; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -t8, 0, t9, t25 * t35 - t5 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41 * t18 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t43, -t45, t24 * t16, t14 - t49, t13, t15, 0, t24 * t36 + t26 * t43, -t24 * t43 + t26 * t36, -t24 * t33 - t5 * t26, t35, -t5 * t24 + t26 * t33, pkin(8) * t35 + t5 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, 0, 0, 0, 0, t16, -t46, t16, t37, t46, pkin(8) * t37 - t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t17, 0.2e1 * t47, 0, 0, 0, pkin(4) * t54, pkin(4) * t55, -0.2e1 * t11 * t26, 0.2e1 * t41 * pkin(8), t11 * t55, t41 * pkin(8) ^ 2 + t11 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t46, t25, t3, -t4, t7 + (t57 - t48) * t25, t32 * t27, 0.2e1 * t39 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t26, t24, -t26, 0, -t24, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t15, -t13, 0, t15, t31 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, 0, -t51, -t50, -t51, t31, t50, t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t57, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
