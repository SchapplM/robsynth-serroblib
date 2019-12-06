% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t34 = cos(qJ(4));
t26 = -t34 * pkin(4) - pkin(3);
t35 = cos(qJ(3));
t41 = t35 * pkin(2);
t18 = t26 - t41;
t48 = 0.2e1 * t18;
t47 = 0.2e1 * t26;
t46 = 0.2e1 * t34;
t29 = sin(qJ(5));
t45 = t29 * pkin(4);
t31 = sin(qJ(3));
t44 = t31 * pkin(2);
t33 = cos(qJ(5));
t43 = t33 * pkin(4);
t42 = t34 * pkin(7);
t25 = -pkin(3) - t41;
t40 = pkin(3) - t25;
t32 = sin(qJ(2));
t36 = cos(qJ(2));
t15 = t31 * t32 - t35 * t36;
t39 = t15 * t34;
t24 = pkin(7) + t44;
t38 = t34 * t24;
t37 = t18 + t26;
t30 = sin(qJ(4));
t16 = t29 * t34 + t33 * t30;
t14 = t29 * t30 - t33 * t34;
t28 = t30 ^ 2;
t27 = t34 * pkin(8);
t21 = t30 * t46;
t20 = t27 + t42;
t19 = (-pkin(7) - pkin(8)) * t30;
t17 = t31 * t36 + t35 * t32;
t13 = t16 ^ 2;
t12 = t27 + t38;
t11 = (-pkin(8) - t24) * t30;
t10 = t15 * t30;
t9 = -t29 * t19 - t33 * t20;
t8 = t33 * t19 - t29 * t20;
t7 = t15 * t16;
t6 = t15 * t14;
t5 = -0.2e1 * t16 * t14;
t4 = -t29 * t11 - t33 * t12;
t3 = t33 * t11 - t29 * t12;
t2 = t14 * t17;
t1 = t16 * t17;
t22 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t36, -t32, 0, -t15, -t17, 0, 0, 0, 0, 0, -t39, t10, 0, 0, 0, 0, 0, t6, t7; 0, 1, 0, 0, 1, 0.2e1 * t41, -0.2e1 * t44, t28, t21, 0, 0, 0, -0.2e1 * t25 * t34, 0.2e1 * t25 * t30, t13, t5, 0, 0, 0, t14 * t48, t16 * t48; 0, 0, 0, 0, 0, -t15, -t17, 0, 0, 0, 0, 0, -t39, t10, 0, 0, 0, 0, 0, t6, t7; 0, 0, 0, 0, 1, t41, -t44, t28, t21, 0, 0, 0, t40 * t34, -t40 * t30, t13, t5, 0, 0, 0, t37 * t14, t37 * t16; 0, 0, 0, 0, 1, 0, 0, t28, t21, 0, 0, 0, pkin(3) * t46, -0.2e1 * pkin(3) * t30, t13, t5, 0, 0, 0, t14 * t47, t16 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30 * t17, -t34 * t17, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, -t30 * t24, -t38, 0, 0, t16, -t14, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, t34, 0, -t30 * pkin(7), -t42, 0, 0, t16, -t14, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t8, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t43, -t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t22;
