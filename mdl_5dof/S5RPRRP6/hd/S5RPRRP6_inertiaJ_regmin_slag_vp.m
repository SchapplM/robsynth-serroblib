% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP6
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
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = sin(qJ(3));
t41 = 0.2e1 * t22;
t23 = cos(qJ(4));
t40 = pkin(3) * t23;
t21 = sin(qJ(4));
t39 = t21 * pkin(4);
t19 = sin(pkin(8));
t10 = t19 * pkin(1) + pkin(6);
t38 = t10 * t21;
t15 = t21 ^ 2;
t37 = t15 * t22;
t36 = t21 * t22;
t35 = t21 * t23;
t24 = cos(qJ(3));
t34 = t21 * t24;
t33 = t23 * t22;
t32 = t23 * t24;
t31 = t24 * t10;
t30 = -qJ(5) - pkin(7);
t29 = qJ(5) * t22;
t28 = t24 * t41;
t27 = t23 * t31;
t20 = cos(pkin(8));
t11 = -t20 * pkin(1) - pkin(2);
t8 = t30 * t21;
t9 = t30 * t23;
t26 = -t8 * t21 - t9 * t23;
t18 = t24 ^ 2;
t17 = t23 ^ 2;
t16 = t22 ^ 2;
t14 = -t23 * pkin(4) - pkin(3);
t13 = t17 * t22;
t12 = t17 * t16;
t7 = -t24 * pkin(3) - t22 * pkin(7) + t11;
t6 = (t10 + t39) * t22;
t5 = t23 * t7;
t4 = t21 * t7 + t27;
t3 = -t21 * t31 + t5;
t2 = t27 + (t7 - t29) * t21;
t1 = -t23 * t29 + t5 + (-pkin(4) - t38) * t24;
t25 = [1, 0, 0, (t19 ^ 2 + t20 ^ 2) * pkin(1) ^ 2, t16, t28, 0, 0, 0, -0.2e1 * t11 * t24, t11 * t41, t12, -0.2e1 * t16 * t35, -0.2e1 * t22 * t32, t21 * t28, t18, 0.2e1 * t16 * t38 - 0.2e1 * t3 * t24, 0.2e1 * t16 * t10 * t23 + 0.2e1 * t4 * t24, (-t1 * t23 - t2 * t21) * t41, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t24 + (-t1 * t21 + t2 * t23) * t22; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15 * t16 + t12 + t18; 0, 0, 0, 0, 0, 0, t22, t24, 0, -t22 * t10, -t31, t21 * t33, t13 - t37, -t34, -t32, 0, -t10 * t33 + (-pkin(3) * t22 + pkin(7) * t24) * t21, pkin(7) * t32 + (t38 - t40) * t22, (-t22 * t8 + t2) * t23 + (t22 * t9 - t1) * t21, t1 * t8 + t6 * t14 - t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0, 0, 0, 0, t32, -t34, t13 + t37, -t24 * t14 + t26 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t15, 0.2e1 * t35, 0, 0, 0, 0.2e1 * t40, -0.2e1 * pkin(3) * t21, 0.2e1 * t26, t14 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t36, -t24, t3, -t4, -pkin(4) * t33, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, -t33, 0, -pkin(4) * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, -t21 * pkin(7), -t23 * pkin(7), -t39, t8 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t25;
