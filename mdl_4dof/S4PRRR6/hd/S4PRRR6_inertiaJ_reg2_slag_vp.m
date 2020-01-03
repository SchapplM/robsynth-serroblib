% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR6_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR6_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR6_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t21 = cos(qJ(3));
t12 = -t21 * pkin(3) - pkin(2);
t33 = 0.2e1 * t12;
t32 = 0.2e1 * t21;
t31 = -pkin(6) - pkin(5);
t17 = sin(qJ(4));
t30 = t17 * pkin(3);
t20 = cos(qJ(4));
t29 = t20 * pkin(3);
t18 = sin(qJ(3));
t19 = sin(qJ(2));
t28 = t18 * t19;
t27 = t21 * t19;
t13 = t18 ^ 2;
t15 = t21 ^ 2;
t26 = t13 + t15;
t25 = t26 * t19;
t7 = t17 * t21 + t20 * t18;
t22 = cos(qJ(2));
t16 = t22 ^ 2;
t14 = t19 ^ 2;
t9 = t31 * t21;
t8 = t31 * t18;
t5 = t17 * t18 - t20 * t21;
t4 = -t17 * t28 + t20 * t27;
t3 = t7 * t19;
t2 = t17 * t8 - t20 * t9;
t1 = t17 * t9 + t20 * t8;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 + t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t14 + t16, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3 ^ 2 + t4 ^ 2 + t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t19, 0, 0, 0, 0, 0, 0, 0, 0, t22 * t21, -t22 * t18, t25, t22 * pkin(2) + pkin(5) * t25, 0, 0, 0, 0, 0, 0, -t22 * t5, -t22 * t7, t3 * t7 - t4 * t5, -t3 * t1 - t22 * t12 + t4 * t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t13, t18 * t32, 0, t15, 0, 0, pkin(2) * t32, -0.2e1 * pkin(2) * t18, 0.2e1 * t26 * pkin(5), t26 * pkin(5) ^ 2 + pkin(2) ^ 2, t7 ^ 2, -0.2e1 * t7 * t5, 0, t5 ^ 2, 0, 0, t5 * t33, t7 * t33, -0.2e1 * t1 * t7 - 0.2e1 * t2 * t5, t1 ^ 2 + t12 ^ 2 + t2 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t27, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, (t17 * t4 - t20 * t3) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t21, 0, -t18 * pkin(5), -t21 * pkin(5), 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, (-t17 * t5 - t20 * t7) * pkin(3), (t1 * t20 + t17 * t2) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t30, 0, (t17 ^ 2 + t20 ^ 2) * pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, -t5, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t30, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
