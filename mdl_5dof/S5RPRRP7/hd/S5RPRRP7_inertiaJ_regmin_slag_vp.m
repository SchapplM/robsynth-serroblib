% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t31 = -t26 * pkin(4) - t24 * qJ(5);
t10 = -pkin(3) + t31;
t49 = -0.2e1 * t10;
t25 = sin(qJ(3));
t48 = 0.2e1 * t25;
t22 = sin(pkin(8));
t12 = t22 * pkin(1) + pkin(6);
t27 = cos(qJ(3));
t37 = t27 * t12;
t23 = cos(pkin(8));
t13 = -t23 * pkin(1) - pkin(2);
t8 = -t27 * pkin(3) - t25 * pkin(7) + t13;
t4 = t24 * t8 + t26 * t37;
t47 = pkin(3) * t24;
t46 = pkin(3) * t26;
t45 = t24 * pkin(7);
t44 = t26 * pkin(7);
t43 = t12 * t24;
t42 = t12 * t26;
t18 = t24 ^ 2;
t41 = t18 * t25;
t40 = t24 * t25;
t39 = t24 * t26;
t38 = t24 * t27;
t16 = t26 * t25;
t17 = t26 * t27;
t20 = t26 ^ 2;
t36 = t18 + t20;
t35 = t27 * qJ(5);
t34 = t27 * t48;
t33 = t36 * pkin(7);
t1 = -t35 + t4;
t7 = t26 * t8;
t2 = -t7 + (pkin(4) + t43) * t27;
t32 = t1 * t26 + t2 * t24;
t30 = -pkin(4) * t24 + t26 * qJ(5);
t21 = t27 ^ 2;
t19 = t25 ^ 2;
t15 = t20 * t25;
t14 = t20 * t19;
t11 = pkin(7) * t38;
t5 = (t12 - t30) * t25;
t3 = -t24 * t37 + t7;
t6 = [1, 0, 0, (t22 ^ 2 + t23 ^ 2) * pkin(1) ^ 2, t19, t34, 0, 0, 0, -0.2e1 * t13 * t27, t13 * t48, t14, -0.2e1 * t19 * t39, -0.2e1 * t25 * t17, t24 * t34, t21, 0.2e1 * t19 * t43 - 0.2e1 * t3 * t27, 0.2e1 * t19 * t42 + 0.2e1 * t4 * t27, 0.2e1 * t2 * t27 + 0.2e1 * t5 * t40, (-t1 * t24 + t2 * t26) * t48, -0.2e1 * t1 * t27 - 0.2e1 * t5 * t16, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32 * t25 - t5 * t27; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18 * t19 + t14 + t21; 0, 0, 0, 0, 0, 0, t25, t27, 0, -t25 * t12, -t37, t24 * t16, t15 - t41, -t38, -t17, 0, t11 + (-t42 - t47) * t25, pkin(7) * t17 + (t43 - t46) * t25, t10 * t40 - t5 * t26 + t11, t32, -t5 * t24 + (-pkin(7) * t27 - t10 * t25) * t26, t32 * pkin(7) + t5 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, 0, 0, 0, 0, t17, -t38, t17, t15 + t41, t38, -t27 * t10 + t25 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t18, 0.2e1 * t39, 0, 0, 0, 0.2e1 * t46, -0.2e1 * t47, t26 * t49, 0.2e1 * t33, t24 * t49, t36 * pkin(7) ^ 2 + t10 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t40, -t27, t3, -t4, t7 + (-0.2e1 * pkin(4) - t43) * t27, t31 * t25, -0.2e1 * t35 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t16, -t40, 0, t16, t30 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, t26, 0, -t45, -t44, -t45, t30, t44, t30 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t16, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
