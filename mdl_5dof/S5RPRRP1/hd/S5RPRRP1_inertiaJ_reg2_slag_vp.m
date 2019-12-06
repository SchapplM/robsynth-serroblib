% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP1
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t34 = cos(qJ(4));
t35 = cos(qJ(3));
t14 = t32 * t35 + t34 * t33;
t12 = t14 ^ 2;
t16 = -t32 * t33 + t34 * t35;
t13 = t16 ^ 2;
t6 = t13 + t12;
t28 = t34 * pkin(3);
t45 = t32 * pkin(3);
t43 = t14 * t45;
t51 = t16 * t28 + t43;
t23 = t28 + pkin(4);
t50 = t16 * t23 + t43;
t21 = t33 * pkin(3) + qJ(2);
t8 = t14 * pkin(4) + t21;
t49 = 0.2e1 * t8;
t48 = 0.2e1 * t21;
t47 = 0.2e1 * qJ(2);
t46 = t16 * pkin(4);
t29 = t33 ^ 2;
t30 = t35 ^ 2;
t20 = t29 + t30;
t36 = -pkin(1) - pkin(6);
t18 = (-pkin(7) + t36) * t33;
t24 = t35 * t36;
t19 = -t35 * pkin(7) + t24;
t4 = -t32 * t18 + t34 * t19;
t2 = -t16 * qJ(5) + t4;
t5 = t34 * t18 + t32 * t19;
t3 = -t14 * qJ(5) + t5;
t41 = t3 * t14 + t2 * t16;
t40 = t5 * t14 + t4 * t16;
t39 = pkin(3) ^ 2;
t38 = 0.2e1 * pkin(4);
t37 = qJ(2) ^ 2;
t26 = t32 ^ 2 * t39;
t25 = -0.2e1 * t45;
t17 = t20 * t36;
t7 = -0.2e1 * t16 * t14;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, -2 * pkin(1), t47, (pkin(1) ^ 2) + t37, t30, -0.2e1 * t35 * t33, 0, t29, 0, 0, t33 * t47, t35 * t47, -0.2e1 * t17, t20 * t36 ^ 2 + t37, t13, t7, 0, t12, 0, 0, t14 * t48, t16 * t48, -0.2e1 * t40, t21 ^ 2 + t4 ^ 2 + t5 ^ 2, t13, t7, 0, t12, 0, 0, t14 * t49, t16 * t49, -0.2e1 * t41, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, -t20, t17, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t40, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, -t33, 0, t24, -t33 * t36, 0, 0, 0, 0, t16, 0, -t14, 0, t4, -t5, -t51, (t32 * t5 + t34 * t4) * pkin(3), 0, 0, t16, 0, -t14, 0, t2, -t3, -t50, t2 * t23 + t3 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t51, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t28, t25, 0, t34 ^ 2 * t39 + t26, 0, 0, 0, 0, 0, 1, 0.2e1 * t23, t25, 0, t23 ^ 2 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, t4, -t5, 0, 0, 0, 0, t16, 0, -t14, 0, t2, -t3, -t46, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t14, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t28, -t45, 0, 0, 0, 0, 0, 0, 0, 1, t38 + t28, -t45, 0, t23 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t38, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
