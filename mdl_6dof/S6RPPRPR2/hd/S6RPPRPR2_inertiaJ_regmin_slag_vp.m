% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
% 
% Output:
% MM_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t31 = cos(pkin(9));
t23 = -t31 * pkin(1) - pkin(2);
t30 = cos(pkin(10));
t19 = -t30 * pkin(3) + t23;
t28 = sin(pkin(10));
t33 = sin(qJ(4));
t47 = cos(qJ(4));
t17 = t47 * t28 + t33 * t30;
t40 = t17 * qJ(5);
t37 = t19 - t40;
t16 = t33 * t28 - t47 * t30;
t49 = t16 * pkin(4);
t6 = t37 + t49;
t53 = -0.2e1 * t6;
t15 = t17 ^ 2;
t52 = 0.2e1 * t19;
t51 = 0.2e1 * qJ(5);
t50 = pkin(4) + pkin(8);
t29 = sin(pkin(9));
t21 = t29 * pkin(1) + qJ(3);
t48 = pkin(7) + t21;
t46 = t16 * t17;
t32 = sin(qJ(6));
t45 = t32 * t16;
t44 = t32 * t17;
t34 = cos(qJ(6));
t10 = t34 * t16;
t11 = t34 * t17;
t43 = t34 * t32;
t42 = t28 ^ 2 + t30 ^ 2;
t41 = qJ(5) * t16;
t39 = 0.2e1 * t46;
t38 = t17 * t50 + t41;
t12 = t48 * t28;
t13 = t48 * t30;
t7 = t47 * t12 + t33 * t13;
t8 = -t33 * t12 + t47 * t13;
t27 = t34 ^ 2;
t26 = t32 ^ 2;
t14 = t16 ^ 2;
t5 = -t16 * pkin(5) + t8;
t4 = t17 * pkin(5) + t7;
t3 = t50 * t16 + t37;
t2 = t34 * t3 + t32 * t4;
t1 = -t32 * t3 + t34 * t4;
t9 = [1, 0, 0 (t29 ^ 2 + t31 ^ 2) * pkin(1) ^ 2, -0.2e1 * t23 * t30, 0.2e1 * t23 * t28, 0.2e1 * t42 * t21, t42 * t21 ^ 2 + t23 ^ 2, t15, -0.2e1 * t46, 0, 0, 0, t16 * t52, t17 * t52, -0.2e1 * t8 * t16 + 0.2e1 * t7 * t17, t16 * t53, t17 * t53, t6 ^ 2 + t7 ^ 2 + t8 ^ 2, t26 * t14, 0.2e1 * t14 * t43, t32 * t39, t34 * t39, t15, 0.2e1 * t1 * t17 - 0.2e1 * t5 * t10, -0.2e1 * t2 * t17 + 0.2e1 * t5 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 * t16 + t8 * t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t15, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t30, t28, 0, t23, 0, 0, 0, 0, 0, t16, t17, 0, -t16, -t17, t6, 0, 0, 0, 0, 0, -t44, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, -t7, -t8, -pkin(4) * t17 - t41, t7, t8, -t7 * pkin(4) + t8 * qJ(5), t16 * t43 (-t26 + t27) * t16, t11, -t44, 0, t5 * t32 - t34 * t38, t32 * t38 + t5 * t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17, 0, t16, t17, t40 - t49, 0, 0, 0, 0, 0, t44, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t51, pkin(4) ^ 2 + qJ(5) ^ 2, t27, -0.2e1 * t43, 0, 0, 0, t32 * t51, t34 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, 0, 0, t7, 0, 0, 0, 0, 0, t11, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, t10, t17, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t32, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32, 0, -t34 * t50, t32 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t9;
