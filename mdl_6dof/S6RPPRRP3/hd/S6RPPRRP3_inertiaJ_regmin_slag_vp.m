% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:04
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:51:54
% EndTime: 2019-05-05 14:51:55
% DurationCPUTime: 0.46s
% Computational Cost: add. (296->74), mult. (505->126), div. (0->0), fcn. (478->6), ass. (0->56)
t60 = 2 * pkin(5);
t27 = sin(pkin(9));
t14 = t27 * pkin(1) + qJ(3);
t59 = 0.2e1 * t14;
t29 = sin(qJ(5));
t58 = -0.2e1 * t29;
t31 = cos(qJ(5));
t57 = 0.2e1 * t31;
t32 = cos(qJ(4));
t56 = 0.2e1 * t32;
t28 = cos(pkin(9));
t16 = -t28 * pkin(1) - pkin(2);
t13 = -pkin(7) + t16;
t30 = sin(qJ(4));
t21 = t31 * t30;
t8 = t30 * pkin(4) - t32 * pkin(8) + t14;
t4 = t13 * t21 + t29 * t8;
t55 = pkin(8) * t30;
t54 = t29 * pkin(8);
t53 = t31 * pkin(8);
t52 = t13 * t29;
t23 = t29 ^ 2;
t51 = t23 * t32;
t50 = t29 * t31;
t18 = t29 * t32;
t49 = t30 * t13;
t22 = t31 * t32;
t36 = -t31 * pkin(5) - t29 * qJ(6);
t12 = -pkin(4) + t36;
t48 = t32 * t12;
t47 = t32 * t13;
t46 = t32 * t30;
t25 = t31 ^ 2;
t45 = t23 + t25;
t24 = t30 ^ 2;
t26 = t32 ^ 2;
t44 = t24 + t26;
t43 = t30 * qJ(6);
t42 = -0.2e1 * t46;
t41 = t45 * pkin(8);
t40 = t45 * t30;
t39 = -pkin(4) * t32 - t55;
t1 = t43 + t4;
t7 = t31 * t8;
t2 = -t7 + (-pkin(5) + t52) * t30;
t38 = t1 * t31 + t2 * t29;
t37 = -t48 + t55;
t35 = -pkin(5) * t29 + t31 * qJ(6);
t20 = t25 * t32;
t19 = t25 * t26;
t17 = t29 * t30;
t11 = t44 * t31;
t10 = t44 * t29;
t5 = (-t13 - t35) * t32;
t3 = -t29 * t49 + t7;
t6 = [1, 0, 0 (t27 ^ 2 + t28 ^ 2) * pkin(1) ^ 2, 0.2e1 * t16, t59, t14 ^ 2 + t16 ^ 2, t26, t42, 0, 0, 0, t30 * t59, t14 * t56, t19, -0.2e1 * t26 * t50, t46 * t57, t29 * t42, t24, -0.2e1 * t26 * t52 + 0.2e1 * t3 * t30, -0.2e1 * t26 * t13 * t31 - 0.2e1 * t4 * t30, 0.2e1 * t5 * t18 - 0.2e1 * t2 * t30 (-t1 * t29 + t2 * t31) * t56, 0.2e1 * t1 * t30 - 0.2e1 * t5 * t22, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5 * t30 + t32 * t38; 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t26 + t19 + t24; 0, 0, 0, 0, 1, 0, t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, -t10, 0, t11, t30 * t38 - t5 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (-0.1e1 + t45) * t46; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45 * t24 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, t47, -t49, t29 * t22, t20 - t51, t17, t21, 0, t39 * t29 + t31 * t47, -t29 * t47 + t31 * t39, -t29 * t37 - t5 * t31, t38, -t5 * t29 + t31 * t37, pkin(8) * t38 + t5 * t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t32, 0, 0, 0, 0, 0, -t21, t17, -t21, t20 + t51, -t17, t30 * t12 + t32 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, 0, 0, 0, 0, t22, -t18, t22, t40, t18, pkin(8) * t40 - t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0.2e1 * t50, 0, 0, 0, pkin(4) * t57, pkin(4) * t58, -0.2e1 * t12 * t31, 0.2e1 * t41, t12 * t58, t45 * pkin(8) ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t18, t30, t3, -t4, t7 + (t60 - t52) * t30, t36 * t32, 0.2e1 * t43 + t4, -t2 * pkin(5) + t1 * qJ(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18, -t22, -t18, 0, t22, t35 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t17, -t21, -t17, 0, t21, t35 * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t31, 0, -t54, -t53, -t54, t35, t53, t35 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t60, 0, 0.2e1 * qJ(6) (pkin(5) ^ 2) + qJ(6) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t22, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t6;
