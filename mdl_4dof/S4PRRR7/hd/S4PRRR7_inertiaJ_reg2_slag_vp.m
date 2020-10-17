% Calculate inertial parameters regressor of joint inertia matrix for
% S4PRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((4+1)*4/2)x(4*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4PRRR7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR7_inertiaJ_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4PRRR7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:36:48
% EndTime: 2019-12-31 16:36:49
% DurationCPUTime: 0.28s
% Computational Cost: add. (133->56), mult. (369->115), div. (0->0), fcn. (400->8), ass. (0->44)
t18 = cos(pkin(4));
t20 = sin(qJ(3));
t23 = cos(qJ(3));
t17 = sin(pkin(4));
t21 = sin(qJ(2));
t39 = t17 * t21;
t5 = -t18 * t23 + t20 * t39;
t46 = t5 ^ 2;
t45 = -0.2e1 * t20;
t44 = 0.2e1 * t23;
t22 = cos(qJ(4));
t43 = pkin(3) * t22;
t14 = t20 ^ 2;
t42 = t14 * pkin(6);
t41 = t20 * pkin(6);
t40 = t5 * t20;
t24 = cos(qJ(2));
t38 = t17 * t24;
t19 = sin(qJ(4));
t37 = t19 * t20;
t36 = t19 * t22;
t35 = t19 * t23;
t34 = t22 * t20;
t33 = t22 * t23;
t13 = t19 ^ 2;
t15 = t22 ^ 2;
t32 = t13 + t15;
t31 = t20 * t44;
t30 = t19 * t34;
t7 = t18 * t20 + t23 * t39;
t1 = -t7 * t19 - t22 * t38;
t2 = -t19 * t38 + t7 * t22;
t29 = -t1 * t19 + t2 * t22;
t8 = -t23 * pkin(3) - t20 * pkin(7) - pkin(2);
t3 = -pkin(6) * t35 + t22 * t8;
t4 = pkin(6) * t33 + t19 * t8;
t28 = -t3 * t19 + t4 * t22;
t27 = t7 * t23 + t40;
t26 = pkin(6) ^ 2;
t16 = t23 ^ 2;
t12 = t17 ^ 2;
t11 = t14 * t26;
t9 = t12 * t24 ^ 2;
t6 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12 * t21 ^ 2 + t18 ^ 2 + t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 ^ 2 + t46 + t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t39, 0, 0, 0, 0, 0, 0, 0, 0, t23 * t38, -t20 * t38, t27, pkin(2) * t38 + t27 * pkin(6), 0, 0, 0, 0, 0, 0, -t1 * t23 + t5 * t37, t2 * t23 + t5 * t34, (-t1 * t22 - t19 * t2) * t20, pkin(6) * t40 + t1 * t3 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t14, t31, 0, t16, 0, 0, pkin(2) * t44, pkin(2) * t45, 0.2e1 * (t14 + t16) * pkin(6), pkin(2) ^ 2 + t16 * t26 + t11, t15 * t14, -0.2e1 * t14 * t36, t33 * t45, t13 * t14, t19 * t31, t16, 0.2e1 * t19 * t42 - 0.2e1 * t3 * t23, 0.2e1 * t22 * t42 + 0.2e1 * t4 * t23, 0.2e1 * (-t19 * t4 - t22 * t3) * t20, t3 ^ 2 + t4 ^ 2 + t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t7, 0, 0, 0, 0, 0, 0, 0, 0, -t5 * t22, t5 * t19, t29, -t5 * pkin(3) + t29 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t23, 0, -t41, -t23 * pkin(6), 0, 0, t30, (-t13 + t15) * t20, -t35, -t30, -t33, 0, -pkin(6) * t34 + (-pkin(3) * t20 + pkin(7) * t23) * t19, pkin(7) * t33 + (pkin(6) * t19 - t43) * t20, t28, -pkin(3) * t41 + t28 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t13, 0.2e1 * t36, 0, t15, 0, 0, 0.2e1 * t43, -0.2e1 * pkin(3) * t19, 0.2e1 * t32 * pkin(7), t32 * pkin(7) ^ 2 + pkin(3) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, -t37, -t23, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, 0, t22, 0, -t19 * pkin(7), -t22 * pkin(7), 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t6;
