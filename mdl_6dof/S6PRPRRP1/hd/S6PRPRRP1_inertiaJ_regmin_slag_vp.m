% Calculate minimal parameter regressor of joint inertia matrix for
% S6PRPRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6PRPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:32:30
% EndTime: 2019-05-04 23:32:31
% DurationCPUTime: 0.36s
% Computational Cost: add. (297->82), mult. (675->146), div. (0->0), fcn. (777->10), ass. (0->54)
t32 = sin(qJ(4));
t56 = 0.2e1 * t32;
t34 = cos(qJ(5));
t55 = pkin(4) * t34;
t31 = sin(qJ(5));
t54 = t31 * pkin(5);
t27 = sin(pkin(11));
t18 = t27 * pkin(2) + pkin(8);
t53 = t18 * t31;
t23 = t31 ^ 2;
t52 = t23 * t32;
t28 = sin(pkin(6));
t33 = sin(qJ(2));
t51 = t28 * t33;
t36 = cos(qJ(2));
t50 = t28 * t36;
t49 = t31 * t32;
t48 = t31 * t34;
t35 = cos(qJ(4));
t47 = t31 * t35;
t46 = t34 * t32;
t45 = t34 * t35;
t44 = t35 * t18;
t43 = -qJ(6) - pkin(9);
t42 = qJ(6) * t32;
t41 = t35 * t56;
t40 = t34 * t44;
t29 = cos(pkin(11));
t19 = -t29 * pkin(2) - pkin(3);
t11 = (t27 * t36 + t29 * t33) * t28;
t30 = cos(pkin(6));
t8 = t11 * t35 + t30 * t32;
t9 = t27 * t51 - t29 * t50;
t1 = -t8 * t31 + t9 * t34;
t2 = t9 * t31 + t8 * t34;
t39 = -t1 * t31 + t2 * t34;
t15 = t43 * t31;
t16 = t43 * t34;
t38 = -t15 * t31 - t16 * t34;
t26 = t35 ^ 2;
t25 = t34 ^ 2;
t24 = t32 ^ 2;
t22 = -t34 * pkin(5) - pkin(4);
t21 = t25 * t32;
t20 = t25 * t24;
t14 = -t35 * pkin(4) - t32 * pkin(9) + t19;
t13 = (t18 + t54) * t32;
t12 = t34 * t14;
t7 = t11 * t32 - t30 * t35;
t6 = t31 * t14 + t40;
t5 = -t31 * t44 + t12;
t4 = t40 + (t14 - t42) * t31;
t3 = -t34 * t42 + t12 + (-pkin(5) - t53) * t35;
t10 = [1, 0, 0, 0, t11 ^ 2 + t30 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, t50, -t51 (t11 * t27 - t29 * t9) * pkin(2), 0, 0, 0, 0, 0, -t9 * t35, t9 * t32, 0, 0, 0, 0, 0, -t1 * t35 + t7 * t49, t2 * t35 + t7 * t46 (-t1 * t34 - t2 * t31) * t32, t1 * t3 + t7 * t13 + t2 * t4; 0, 1, 0, 0 (t27 ^ 2 + t29 ^ 2) * pkin(2) ^ 2, t24, t41, 0, 0, 0, -0.2e1 * t19 * t35, t19 * t56, t20, -0.2e1 * t24 * t48, -0.2e1 * t32 * t45, t31 * t41, t26, 0.2e1 * t24 * t53 - 0.2e1 * t5 * t35, 0.2e1 * t24 * t18 * t34 + 0.2e1 * t6 * t35 (-t3 * t34 - t31 * t4) * t56, t13 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, t30, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * t32 - t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13 * t35 + (-t3 * t31 + t4 * t34) * t32; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t24 + t20 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, -t8, 0, 0, 0, 0, 0, -t7 * t34, t7 * t31, t39, t1 * t15 - t2 * t16 + t7 * t22; 0, 0, 0, 0, 0, 0, 0, t32, t35, 0, -t32 * t18, -t44, t31 * t46, t21 - t52, -t47, -t45, 0, -t18 * t46 + (-pkin(4) * t32 + pkin(9) * t35) * t31, pkin(9) * t45 + (t53 - t55) * t32 (-t15 * t32 + t4) * t34 + (t16 * t32 - t3) * t31, t13 * t22 + t3 * t15 - t4 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t32, 0, 0, 0, 0, 0, t45, -t47, t21 + t52, -t35 * t22 + t38 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0.2e1 * t48, 0, 0, 0, 0.2e1 * t55, -0.2e1 * pkin(4) * t31, 0.2e1 * t38, t15 ^ 2 + t16 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t46, -t49, -t35, t5, -t6, -pkin(5) * t46, t3 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t49, -t46, 0, -pkin(5) * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t34, 0, -t31 * pkin(9), -t34 * pkin(9), -t54, t15 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t10;
