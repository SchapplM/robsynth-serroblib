% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d4,d5,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x19]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRPRR4_inertiaJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(4));
t39 = 0.2e1 * t22;
t24 = cos(qJ(5));
t38 = pkin(4) * t24;
t17 = sin(pkin(10));
t12 = t17 * pkin(2) + pkin(7);
t21 = sin(qJ(5));
t37 = t12 * t21;
t18 = sin(pkin(5));
t23 = sin(qJ(2));
t36 = t18 * t23;
t26 = cos(qJ(2));
t35 = t18 * t26;
t34 = t21 * t22;
t33 = t21 * t24;
t25 = cos(qJ(4));
t32 = t21 * t25;
t31 = t24 * t22;
t30 = t24 * t25;
t29 = t25 * t12;
t28 = t25 * t39;
t19 = cos(pkin(10));
t13 = -t19 * pkin(2) - pkin(3);
t20 = cos(pkin(5));
t16 = t24 ^ 2;
t15 = t22 ^ 2;
t14 = t21 ^ 2;
t10 = -t25 * pkin(4) - t22 * pkin(8) + t13;
t9 = (t17 * t26 + t19 * t23) * t18;
t7 = t17 * t36 - t19 * t35;
t6 = t20 * t22 + t9 * t25;
t5 = -t20 * t25 + t9 * t22;
t4 = t21 * t10 + t24 * t29;
t3 = t24 * t10 - t21 * t29;
t2 = t7 * t21 + t6 * t24;
t1 = -t6 * t21 + t7 * t24;
t8 = [1, 0, 0, 0, t20 ^ 2 + t7 ^ 2 + t9 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, t35, -t36, (t17 * t9 - t19 * t7) * pkin(2), 0, 0, 0, 0, 0, -t7 * t25, t7 * t22, 0, 0, 0, 0, 0, -t1 * t25 + t5 * t34, t2 * t25 + t5 * t31; 0, 1, 0, 0, (t17 ^ 2 + t19 ^ 2) * pkin(2) ^ 2, t15, t28, 0, 0, 0, -0.2e1 * t13 * t25, t13 * t39, t16 * t15, -0.2e1 * t15 * t33, -0.2e1 * t22 * t30, t21 * t28, t25 ^ 2, 0.2e1 * t15 * t37 - 0.2e1 * t3 * t25, 0.2e1 * t15 * t12 * t24 + 0.2e1 * t4 * t25; 0, 0, 0, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t5 * t24, t5 * t21; 0, 0, 0, 0, 0, 0, 0, t22, t25, 0, -t22 * t12, -t29, t21 * t31, (-t14 + t16) * t22, -t32, -t30, 0, -t12 * t31 + (-pkin(4) * t22 + pkin(8) * t25) * t21, pkin(8) * t30 + (t37 - t38) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t22, 0, 0, 0, 0, 0, t30, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t14, 0.2e1 * t33, 0, 0, 0, 0.2e1 * t38, -0.2e1 * pkin(4) * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t34, -t25, t3, -t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t34, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t24, 0, -t21 * pkin(8), -t24 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
