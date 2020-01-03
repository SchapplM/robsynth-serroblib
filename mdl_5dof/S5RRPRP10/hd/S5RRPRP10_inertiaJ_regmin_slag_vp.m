% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP10_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(2));
t47 = -0.2e1 * t26;
t28 = cos(qJ(2));
t46 = 0.2e1 * t28;
t45 = 2 * qJ(3);
t29 = -pkin(2) - pkin(7);
t27 = cos(qJ(4));
t44 = t27 * pkin(4);
t18 = t26 * pkin(6);
t12 = t26 * pkin(3) + t18;
t25 = sin(qJ(4));
t43 = t25 * t12;
t42 = t25 * t26;
t41 = t25 * t28;
t40 = t26 * t28;
t39 = t27 * t25;
t38 = t27 * t28;
t19 = t28 * pkin(6);
t13 = t28 * pkin(3) + t19;
t22 = t26 ^ 2;
t24 = t28 ^ 2;
t37 = t22 + t24;
t36 = t28 * qJ(3);
t35 = -0.2e1 * t40;
t34 = -t26 * qJ(3) - pkin(1);
t7 = t29 * t28 + t34;
t33 = qJ(5) * t28 - t7;
t32 = -t26 * pkin(2) + t36;
t31 = t26 * t29 + t36;
t23 = t27 ^ 2;
t21 = t25 ^ 2;
t17 = t27 * t29;
t16 = t27 * t26;
t15 = t25 * pkin(4) + qJ(3);
t14 = t21 + t23;
t11 = -t28 * pkin(2) + t34;
t10 = -t27 * qJ(5) + t17;
t9 = (-qJ(5) + t29) * t25;
t8 = t27 * t12;
t6 = pkin(4) * t38 + t13;
t5 = t10 * t27 + t9 * t25;
t4 = t27 * t7 + t43;
t3 = -t25 * t7 + t8;
t2 = -t33 * t27 + t43;
t1 = t26 * pkin(4) + t33 * t25 + t8;
t20 = [1, 0, 0, t22, 0.2e1 * t40, 0, 0, 0, pkin(1) * t46, pkin(1) * t47, 0.2e1 * t37 * pkin(6), t11 * t46, t11 * t47, t37 * pkin(6) ^ 2 + t11 ^ 2, t21 * t24, 0.2e1 * t24 * t39, t25 * t35, t27 * t35, t22, 0.2e1 * t13 * t38 + 0.2e1 * t3 * t26, -0.2e1 * t13 * t41 - 0.2e1 * t4 * t26, (t1 * t25 - t2 * t27) * t46, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t26, t28, 0, -t18, -t19, t32, t18, t19, t32 * pkin(6), -t25 * t38, (t21 - t23) * t28, t16, -t42, 0, t13 * t25 + t31 * t27, t13 * t27 - t31 * t25, (-t28 * t9 - t1) * t27 + (t10 * t28 - t2) * t25, t1 * t10 + t6 * t15 + t2 * t9; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t45, pkin(2) ^ 2 + (qJ(3) ^ 2), t23, -0.2e1 * t39, 0, 0, 0, t25 * t45, t27 * t45, -0.2e1 * t5, t10 ^ 2 + t15 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, 0, t18, 0, 0, 0, 0, 0, t16, -t42, 0, t1 * t27 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, -t14, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t38, t26, t3, -t4, pkin(4) * t41, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t17, -t25 * t29, -t44, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t25, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t20;
