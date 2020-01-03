% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(3));
t29 = cos(qJ(3));
t30 = -pkin(1) - pkin(2);
t15 = t26 * qJ(2) - t29 * t30;
t13 = pkin(3) + t15;
t28 = cos(qJ(4));
t41 = t28 * pkin(4);
t9 = t13 + t41;
t49 = -0.2e1 * t9;
t24 = sin(qJ(5));
t25 = sin(qJ(4));
t27 = cos(qJ(5));
t11 = t24 * t28 + t27 * t25;
t48 = t11 ^ 2;
t21 = -pkin(3) - t41;
t47 = 0.2e1 * t21;
t46 = -0.2e1 * t25;
t45 = 0.2e1 * t28;
t44 = pkin(7) + pkin(8);
t43 = t24 * pkin(4);
t42 = t27 * pkin(4);
t40 = pkin(3) + t13;
t16 = t29 * qJ(2) + t26 * t30;
t14 = -pkin(7) + t16;
t39 = pkin(8) - t14;
t38 = -t21 + t9;
t10 = t24 * t25 - t27 * t28;
t37 = t11 * t10;
t36 = t25 * t28;
t35 = t29 * t10;
t34 = t29 * t11;
t33 = t29 * t25;
t32 = t29 * t28;
t31 = -0.2e1 * t37;
t23 = t25 ^ 2;
t19 = 0.2e1 * t36;
t18 = t44 * t28;
t17 = t44 * t25;
t8 = t10 * t26;
t7 = t11 * t26;
t6 = t39 * t28;
t5 = t39 * t25;
t4 = t24 * t17 - t27 * t18;
t3 = -t27 * t17 - t24 * t18;
t2 = -t24 * t5 + t27 * t6;
t1 = t24 * t6 + t27 * t5;
t12 = [1, 0, 0, 2 * pkin(1), 0.2e1 * qJ(2), (pkin(1) ^ 2) + qJ(2) ^ 2, 1, 0.2e1 * t15, 0.2e1 * t16, t23, t19, 0, 0, 0, t13 * t45, t13 * t46, t48, t31, 0, 0, 0, t10 * t49, t11 * t49; 0, 0, 0, -1, 0, -pkin(1), 0, -t29, t26, 0, 0, 0, 0, 0, -t32, t33, 0, 0, 0, 0, 0, t35, t34; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, -1, -t15, -t16, -t23, -0.2e1 * t36, 0, 0, 0, -t40 * t28, t40 * t25, -t48, 0.2e1 * t37, 0, 0, 0, t38 * t10, t38 * t11; 0, 0, 0, 0, 0, 0, 0, t29, -t26, 0, 0, 0, 0, 0, t32, -t33, 0, 0, 0, 0, 0, -t35, -t34; 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, t19, 0, 0, 0, pkin(3) * t45, pkin(3) * t46, t48, t31, 0, 0, 0, t10 * t47, t11 * t47; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t28, 0, -t25 * t14, -t28 * t14, 0, 0, -t11, t10, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25 * t26, -t28 * t26, 0, 0, 0, 0, 0, -t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * pkin(7), -t28 * pkin(7), 0, 0, t11, -t10, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, t10, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, -t10, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t42, -t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
