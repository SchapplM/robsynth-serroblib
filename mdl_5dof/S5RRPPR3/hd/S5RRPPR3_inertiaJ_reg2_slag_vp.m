% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:26
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR3_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t28 = cos(qJ(2));
t20 = t28 * pkin(1);
t19 = t20 + pkin(2);
t23 = sin(pkin(8));
t24 = cos(pkin(8));
t26 = sin(qJ(2));
t38 = t26 * pkin(1);
t32 = t24 * t38;
t45 = -(pkin(2) + t19) * t23 - t32;
t9 = t23 * t19 + t32;
t5 = qJ(4) + t9;
t44 = t5 ^ 2;
t40 = t23 * pkin(2);
t16 = qJ(4) + t40;
t43 = t16 ^ 2;
t25 = sin(qJ(5));
t42 = 0.2e1 * t25;
t27 = cos(qJ(5));
t41 = 0.2e1 * t27;
t39 = t24 * pkin(2);
t37 = t5 * t16;
t35 = t16 + t5;
t34 = -t24 * t19 + t23 * t38;
t21 = t25 ^ 2;
t22 = t27 ^ 2;
t33 = t21 + t22;
t7 = -pkin(3) + t34;
t18 = -pkin(3) - t39;
t14 = -pkin(7) + t18;
t2 = t33 * t14;
t31 = -t34 + t39;
t15 = -0.2e1 * t27 * t25;
t4 = -pkin(7) + t7;
t1 = t33 * t4;
t3 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t38, 0, (t26 ^ 2 + t28 ^ 2) * pkin(1) ^ 2, 0, 0, 0, 0, 0, 1, -0.2e1 * t34, -0.2e1 * t9, 0, t34 ^ 2 + t9 ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t7, 0.2e1 * t5, t7 ^ 2 + t44, t22, t15, 0, t21, 0, 0, t5 * t42, t5 * t41, -0.2e1 * t1, t33 * t4 ^ 2 + t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t20, -t38, 0, 0, 0, 0, 0, 0, 0, 1, t31, t45, 0, (t23 * t9 - t24 * t34) * pkin(2), 1, 0, 0, 0, 0, 0, 0, -0.2e1 * pkin(3) - t31, 0.2e1 * qJ(4) - t45, t7 * t18 + t37, t22, t15, 0, t21, 0, 0, t35 * t25, t35 * t27, (-t14 - t4) * t33, t2 * t4 + t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t39, -0.2e1 * t40, 0, (t23 ^ 2 + t24 ^ 2) * pkin(2) ^ 2, 1, 0, 0, 0, 0, 0, 0, 0.2e1 * t18, 0.2e1 * t16, t18 ^ 2 + t43, t22, t15, 0, t21, 0, 0, t16 * t42, t16 * t41, -0.2e1 * t2, t14 ^ 2 * t33 + t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t7, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t18, 0, 0, 0, 0, 0, 0, 0, 0, -t33, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t27 * t4, -t25 * t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, 0, -t25, 0, t27 * t14, -t25 * t14, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t27, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t3;
