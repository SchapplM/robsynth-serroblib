% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP13_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP13_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t50 = 2 * pkin(4);
t20 = sin(qJ(4));
t49 = -0.2e1 * t20;
t22 = cos(qJ(4));
t48 = 0.2e1 * t22;
t47 = 2 * qJ(2);
t21 = sin(qJ(3));
t46 = pkin(7) * t21;
t45 = t20 * pkin(7);
t44 = t22 * pkin(7);
t23 = cos(qJ(3));
t10 = t21 * pkin(3) - t23 * pkin(7) + qJ(2);
t24 = -pkin(1) - pkin(6);
t39 = t22 * t24;
t4 = t20 * t10 + t21 * t39;
t43 = t20 * t22;
t42 = t20 * t23;
t41 = t20 * t24;
t40 = t21 * t24;
t15 = t22 * t23;
t27 = -t22 * pkin(4) - t20 * qJ(5);
t11 = -pkin(3) + t27;
t38 = t23 * t11;
t37 = t23 * t21;
t36 = t23 * t24;
t16 = t20 ^ 2;
t18 = t22 ^ 2;
t35 = t16 + t18;
t17 = t21 ^ 2;
t19 = t23 ^ 2;
t34 = t17 + t19;
t33 = t21 * qJ(5);
t32 = -0.2e1 * t37;
t31 = t35 * t21;
t30 = -pkin(3) * t23 - t46;
t1 = t33 + t4;
t7 = t22 * t10;
t2 = -t7 + (-pkin(4) + t41) * t21;
t29 = t1 * t22 + t2 * t20;
t28 = -t38 + t46;
t26 = -pkin(4) * t20 + t22 * qJ(5);
t14 = t22 * t21;
t13 = t20 * t21;
t9 = t34 * t22;
t8 = t34 * t20;
t5 = (-t24 - t26) * t23;
t3 = -t20 * t40 + t7;
t6 = [1, 0, 0, -2 * pkin(1), t47, pkin(1) ^ 2 + qJ(2) ^ 2, t19, t32, 0, 0, 0, t21 * t47, t23 * t47, t18 * t19, -0.2e1 * t19 * t43, t37 * t48, t20 * t32, t17, -0.2e1 * t19 * t41 + 0.2e1 * t3 * t21, -0.2e1 * t19 * t39 - 0.2e1 * t4 * t21, -0.2e1 * t2 * t21 + 0.2e1 * t5 * t42, 0.2e1 * (-t1 * t20 + t2 * t22) * t23, 0.2e1 * t1 * t21 - 0.2e1 * t5 * t15, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, -t9, -t8, 0, t9, t29 * t21 - t5 * t23; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35 * t17 + t19; 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, t36, -t40, t20 * t15, (-t16 + t18) * t23, t13, t14, 0, t30 * t20 + t22 * t36, -t20 * t36 + t30 * t22, -t28 * t20 - t5 * t22, t29, -t5 * t20 + t28 * t22, t29 * pkin(7) + t5 * t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, -t21, 0, 0, 0, 0, 0, t15, -t42, t15, t31, t42, pkin(7) * t31 - t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t16, 0.2e1 * t43, 0, 0, 0, pkin(3) * t48, pkin(3) * t49, -0.2e1 * t11 * t22, 0.2e1 * t35 * pkin(7), t11 * t49, t35 * pkin(7) ^ 2 + t11 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t42, t21, t3, -t4, t7 + (t50 - t41) * t21, t27 * t23, 0.2e1 * t33 + t4, -t2 * pkin(4) + t1 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, -t14, -t13, 0, t14, t26 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t22, 0, -t45, -t44, -t45, t26, t44, t26 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t50, 0, 0.2e1 * qJ(5), (pkin(4) ^ 2) + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21, t15, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
