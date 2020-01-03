% Calculate inertial parameters regressor of joint inertia matrix for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRRR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t25 = sin(qJ(4));
t23 = t25 ^ 2;
t28 = cos(qJ(4));
t24 = t28 ^ 2;
t37 = t23 + t24;
t30 = cos(qJ(2));
t21 = t30 * pkin(1);
t17 = t21 + pkin(2);
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t27 = sin(qJ(2));
t42 = t27 * pkin(1);
t36 = t29 * t42;
t8 = t26 * t17 + t36;
t6 = pkin(7) + t8;
t45 = t37 * t6;
t43 = t26 * pkin(2);
t15 = pkin(7) + t43;
t46 = t37 * t15;
t44 = pkin(3) * t25;
t34 = -t29 * t17 + t26 * t42;
t5 = -pkin(3) + t34;
t41 = t5 * t28;
t20 = t29 * pkin(2);
t16 = -t20 - pkin(3);
t39 = t16 * t28;
t38 = pkin(7) * t37;
t22 = pkin(3) * t28;
t13 = 0.2e1 * t25 * t28;
t12 = t16 * t25;
t3 = t5 * t25;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t21, -0.2e1 * t42, 0, (t27 ^ 2 + t30 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t34, -0.2e1 * t8, 0, t34 ^ 2 + t8 ^ 2, t23, t13, 0, t24, 0, 0, -0.2e1 * t41, 0.2e1 * t3, 0.2e1 * t45, t37 * t6 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t21, -t42, 0, 0, 0, 0, 0, 0, 0, 1, t20 - t34, -t36 + (-pkin(2) - t17) * t26, 0, (t26 * t8 - t29 * t34) * pkin(2), t23, t13, 0, t24, 0, 0, (-t16 - t5) * t28, t12 + t3, t46 + t45, t5 * t16 + t46 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t43, 0, (t26 ^ 2 + t29 ^ 2) * pkin(2) ^ 2, t23, t13, 0, t24, 0, 0, -0.2e1 * t39, 0.2e1 * t12, 0.2e1 * t46, t15 ^ 2 * t37 + t16 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, -t34, -t8, 0, 0, t23, t13, 0, t24, 0, 0, t22 - t41, t3 - t44, t38 + t45, -t5 * pkin(3) + pkin(7) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t20, -t43, 0, 0, t23, t13, 0, t24, 0, 0, t22 - t39, t12 - t44, t38 + t46, -t16 * pkin(3) + pkin(7) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t23, t13, 0, t24, 0, 0, 0.2e1 * t22, -0.2e1 * t44, 0.2e1 * t38, pkin(7) ^ 2 * t37 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * t6, -t28 * t6, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * t15, -t28 * t15, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, 0, t28, 0, -t25 * pkin(7), -t28 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
