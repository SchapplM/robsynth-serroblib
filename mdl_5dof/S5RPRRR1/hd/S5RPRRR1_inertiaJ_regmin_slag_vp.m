% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% 
% Output:
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-07-18 13:26
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(1,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S5RPRRR1_inertiaJ_regmin_slag_vp: pkin has to be [1x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = sin(qJ(5));
t14 = cos(qJ(5));
t16 = cos(qJ(3));
t28 = t16 * t14;
t13 = sin(qJ(3));
t15 = cos(qJ(4));
t32 = t15 * t13;
t3 = t11 * t32 + t28;
t47 = -0.2e1 * t3;
t12 = sin(qJ(4));
t46 = 0.2e1 * t12;
t45 = 2 * qJ(2);
t44 = t15 * t3;
t29 = t16 * t11;
t31 = t15 * t14;
t4 = t13 * t31 - t29;
t43 = t4 * t11;
t42 = t4 * t15;
t41 = t11 * t12;
t40 = t11 * t13;
t39 = t11 * t14;
t38 = t12 * t13;
t37 = t12 * t15;
t36 = t12 * t16;
t35 = t13 * t14;
t34 = t14 * t12;
t33 = t15 * t11;
t30 = t15 * t16;
t27 = t13 * qJ(2);
t26 = t16 * qJ(2);
t25 = -0.2e1 * t37;
t24 = 0.2e1 * t37;
t23 = 0.2e1 * t13 * t16;
t6 = t12 ^ 2;
t22 = t6 * t40;
t21 = t6 * t35;
t20 = t6 * t26;
t19 = t12 * t32;
t18 = t12 * t26;
t10 = t16 ^ 2;
t7 = t13 ^ 2;
t17 = (t10 + t7) * t45;
t9 = t15 ^ 2;
t8 = t14 ^ 2;
t5 = t11 ^ 2;
t2 = (t15 * t28 + t40) * qJ(2);
t1 = (-t15 * t29 + t35) * qJ(2);
t48 = [1, 0, 0, 0, t45, qJ(2) ^ 2, t7, t23, 0, 0, 0, 0, 0, t9 * t7, t7 * t25, -0.2e1 * t13 * t30, t12 * t23, t10, t12 * t17, t15 * t17, t4 ^ 2, t4 * t47, 0.2e1 * t4 * t38, t38 * t47, t6 * t7, (t1 * t13 + t3 * t26) * t46, (-t13 * t2 + t4 * t26) * t46; 0, 0, 0, -1, 0, 0, 0, 0, 0, 0, 0, -t16, t13, 0, 0, 0, 0, 0, -t30, t36, 0, 0, 0, 0, 0, -t22 - t44, -t21 - t42; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t13, t16, 0, -t27, -t26, t19, (-t6 + t9) * t13, -t36, -t30, 0, -t15 * t27, t12 * t27, t4 * t34, (-t14 * t3 - t43) * t12, t21 - t42, -t22 + t44, -t19, -t1 * t15 + t11 * t20, t14 * t20 + t2 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t6, t24, 0, 0, 0, 0, 0, t8 * t6, -0.2e1 * t6 * t39, t14 * t25, t11 * t24, t9, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t38, -t16, -t18, -t15 * t26, t43, -t11 * t3 + t4 * t14, t11 * t38, t13 * t34, 0, -t14 * t18, t11 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t12, 0, 0, 0, 0, 0, t31, -t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t15, 0, 0, 0, t11 * t34, (-t5 + t8) * t12, -t33, -t31, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t5, 0.2e1 * t39, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t3, t38, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t41, -t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t14, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t48;
