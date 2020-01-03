% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP10
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
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t25 = sin(pkin(8));
t26 = cos(pkin(8));
t28 = sin(qJ(3));
t44 = cos(qJ(3));
t13 = t44 * t25 + t28 * t26;
t48 = -0.2e1 * t13;
t19 = -t26 * pkin(2) - pkin(1);
t47 = 0.2e1 * t19;
t29 = cos(qJ(4));
t46 = t29 * pkin(4);
t39 = pkin(6) + qJ(2);
t14 = t39 * t25;
t15 = t39 * t26;
t9 = -t28 * t14 + t44 * t15;
t45 = t29 * t9;
t12 = t28 * t25 - t44 * t26;
t27 = sin(qJ(4));
t43 = t27 * t12;
t42 = t27 * t13;
t41 = t27 * t29;
t40 = t29 * t13;
t38 = -qJ(5) - pkin(7);
t37 = t25 ^ 2 + t26 ^ 2;
t23 = t27 ^ 2;
t24 = t29 ^ 2;
t36 = t23 + t24;
t35 = qJ(5) * t13;
t34 = t12 * t48;
t7 = t12 * pkin(3) - t13 * pkin(7) + t19;
t3 = -t27 * t9 + t29 * t7;
t33 = -pkin(3) * t13 - pkin(7) * t12;
t1 = t12 * pkin(4) - t29 * t35 + t3;
t2 = t45 + (t7 - t35) * t27;
t32 = t1 * t29 + t2 * t27;
t16 = t38 * t27;
t17 = t38 * t29;
t31 = t29 * t16 - t27 * t17;
t8 = t44 * t14 + t28 * t15;
t20 = -pkin(3) - t46;
t11 = t13 ^ 2;
t10 = t29 * t12;
t5 = pkin(4) * t42 + t8;
t4 = t27 * t7 + t45;
t6 = [1, 0, 0, 0.2e1 * pkin(1) * t26, -0.2e1 * pkin(1) * t25, 0.2e1 * t37 * qJ(2), t37 * qJ(2) ^ 2 + pkin(1) ^ 2, t11, t34, 0, 0, 0, t12 * t47, t13 * t47, t24 * t11, -0.2e1 * t11 * t41, 0.2e1 * t12 * t40, t27 * t34, t12 ^ 2, 0.2e1 * t3 * t12 + 0.2e1 * t8 * t42, -0.2e1 * t4 * t12 + 0.2e1 * t8 * t40, t32 * t48, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, -t26, t25, 0, -pkin(1), 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, t10, -t43, -t36 * t13, t32; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t12, 0, -t8, -t9, t27 * t40, (-t23 + t24) * t13, t43, t10, 0, t33 * t27 - t8 * t29, t8 * t27 + t33 * t29, -t1 * t27 - t31 * t13 + t2 * t29, t1 * t16 - t2 * t17 + t5 * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0.2e1 * t41, 0, 0, 0, 0.2e1 * pkin(3) * t29, -0.2e1 * pkin(3) * t27, -0.2e1 * t16 * t27 - 0.2e1 * t17 * t29, t16 ^ 2 + t17 ^ 2 + t20 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, -t42, t12, t3, -t4, -pkin(4) * t40, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t27, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t29, 0, -t27 * pkin(7), -t29 * pkin(7), -t27 * pkin(4), t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
