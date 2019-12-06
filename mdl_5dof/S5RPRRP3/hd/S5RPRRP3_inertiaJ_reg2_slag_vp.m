% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRRP3
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
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:04
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t31 = cos(pkin(8));
t42 = t31 * pkin(1);
t22 = -pkin(2) - t42;
t35 = cos(qJ(3));
t19 = -t35 * pkin(3) + t22;
t32 = sin(qJ(4));
t33 = sin(qJ(3));
t34 = cos(qJ(4));
t16 = t32 * t33 - t34 * t35;
t44 = t16 * pkin(4);
t8 = t19 + t44;
t47 = 0.2e1 * t8;
t46 = 0.2e1 * t19;
t45 = 0.2e1 * t33;
t30 = sin(pkin(8));
t43 = t30 * pkin(1);
t41 = t32 * pkin(3);
t27 = t34 * pkin(3);
t21 = pkin(6) + t43;
t40 = pkin(7) + t21;
t28 = t33 ^ 2;
t29 = t35 ^ 2;
t39 = t28 + t29;
t12 = t40 * t33;
t13 = t40 * t35;
t4 = -t34 * t12 - t32 * t13;
t5 = -t32 * t12 + t34 * t13;
t37 = pkin(3) ^ 2;
t36 = 0.2e1 * pkin(4);
t26 = t32 ^ 2 * t37;
t25 = -0.2e1 * t41;
t24 = t27 + pkin(4);
t18 = t32 * t35 + t34 * t33;
t15 = t18 ^ 2;
t14 = t16 ^ 2;
t11 = t18 * t41;
t10 = t16 * t41;
t7 = -0.2e1 * t18 * t16;
t6 = t15 + t14;
t3 = -t16 * qJ(5) + t5;
t2 = -t18 * qJ(5) + t4;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t43, 0, (t30 ^ 2 + t31 ^ 2) * pkin(1) ^ 2, t28, t35 * t45, 0, t29, 0, 0, -0.2e1 * t22 * t35, t22 * t45, 0.2e1 * t39 * t21, t39 * t21 ^ 2 + t22 ^ 2, t15, t7, 0, t14, 0, 0, t16 * t46, t18 * t46, -0.2e1 * t5 * t16 - 0.2e1 * t4 * t18, t19 ^ 2 + t4 ^ 2 + t5 ^ 2, t15, t7, 0, t14, 0, 0, t16 * t47, t18 * t47, -0.2e1 * t3 * t16 - 0.2e1 * t2 * t18, t2 ^ 2 + t3 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4 * t16 + t5 * t18, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2 * t16 + t3 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, t35, 0, -t33 * t21, -t35 * t21, 0, 0, 0, 0, t18, 0, -t16, 0, t4, -t5, -t18 * t27 - t10, (t32 * t5 + t34 * t4) * pkin(3), 0, 0, t18, 0, -t16, 0, t2, -t3, -t24 * t18 - t10, t2 * t24 + t3 * t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, -t33, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, -t16 * t27 + t11, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, -t16 * t24 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t27, t25, 0, t34 ^ 2 * t37 + t26, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, t25, 0, t24 ^ 2 + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, t4, -t5, 0, 0, 0, 0, t18, 0, -t16, 0, t2, -t3, -t18 * pkin(4), t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t18, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t27, -t41, 0, 0, 0, 0, 0, 0, 0, 1, t36 + t27, -t41, 0, t24 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
