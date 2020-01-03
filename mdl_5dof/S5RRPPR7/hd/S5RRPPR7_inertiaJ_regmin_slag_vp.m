% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR7
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
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(pkin(8));
t27 = cos(pkin(8));
t29 = sin(qJ(2));
t31 = cos(qJ(2));
t13 = t26 * t29 - t27 * t31;
t14 = t26 * t31 + t27 * t29;
t23 = -t31 * pkin(2) - pkin(1);
t34 = -t14 * qJ(4) + t23;
t6 = t13 * pkin(3) + t34;
t46 = -0.2e1 * t6;
t19 = t26 * pkin(2) + qJ(4);
t45 = 0.2e1 * t19;
t44 = 0.2e1 * t31;
t43 = t19 * t13;
t28 = sin(qJ(5));
t42 = t28 * t13;
t41 = t28 * t14;
t30 = cos(qJ(5));
t40 = t30 * t13;
t11 = t30 * t14;
t39 = t30 * t28;
t38 = -qJ(3) - pkin(6);
t37 = 0.2e1 * t13 * t14;
t16 = t38 * t29;
t17 = t38 * t31;
t7 = -t27 * t16 - t26 * t17;
t9 = t26 * t16 - t27 * t17;
t36 = t7 ^ 2 + t9 ^ 2;
t22 = -t27 * pkin(2) - pkin(3);
t18 = -pkin(7) + t22;
t35 = -t14 * t18 + t43;
t33 = -0.2e1 * t9 * t13 + 0.2e1 * t7 * t14;
t25 = t30 ^ 2;
t24 = t28 ^ 2;
t12 = t13 ^ 2;
t5 = -t13 * pkin(4) + t9;
t4 = t14 * pkin(4) + t7;
t3 = (pkin(3) + pkin(7)) * t13 + t34;
t2 = t28 * t4 + t30 * t3;
t1 = -t28 * t3 + t30 * t4;
t8 = [1, 0, 0, t29 ^ 2, t29 * t44, 0, 0, 0, pkin(1) * t44, -0.2e1 * pkin(1) * t29, t33, t23 ^ 2 + t36, t33, t13 * t46, t14 * t46, t6 ^ 2 + t36, t24 * t12, 0.2e1 * t12 * t39, t28 * t37, t30 * t37, t14 ^ 2, 0.2e1 * t1 * t14 - 0.2e1 * t5 * t40, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t42; 0, 0, 0, 0, 0, t29, t31, 0, -t29 * pkin(6), -t31 * pkin(6), (-t13 * t26 - t14 * t27) * pkin(2), (t26 * t9 - t27 * t7) * pkin(2), t22 * t14 - t43, t7, t9, t9 * t19 + t7 * t22, t13 * t39, (-t24 + t25) * t13, t11, -t41, 0, t5 * t28 - t35 * t30, t35 * t28 + t5 * t30; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t26 ^ 2 + t27 ^ 2) * pkin(2) ^ 2, 0, 0.2e1 * t22, t45, t19 ^ 2 + t22 ^ 2, t25, -0.2e1 * t39, 0, 0, 0, t28 * t45, t30 * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, -t13, -t14, t6, 0, 0, 0, 0, 0, -t41, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, 0, t7, 0, 0, 0, 0, 0, t11, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, t22, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, t40, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28, 0, t30 * t18, -t28 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
