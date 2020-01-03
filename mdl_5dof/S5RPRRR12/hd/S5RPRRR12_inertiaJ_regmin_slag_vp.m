% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR12
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
% MM_reg [((5+1)*5/2)x27]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR12_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t29 = sin(qJ(4));
t30 = sin(qJ(3));
t32 = cos(qJ(3));
t41 = cos(qJ(4));
t16 = t29 * t32 + t41 * t30;
t13 = t16 ^ 2;
t15 = t29 * t30 - t41 * t32;
t14 = t15 ^ 2;
t48 = -t13 - t14;
t47 = -0.2e1 * t15;
t46 = 0.2e1 * t16;
t45 = 2 * qJ(2);
t44 = t29 * pkin(3);
t31 = cos(qJ(5));
t33 = -pkin(1) - pkin(6);
t18 = (-pkin(7) + t33) * t30;
t24 = t32 * t33;
t36 = -t32 * pkin(7) + t24;
t6 = t29 * t18 - t41 * t36;
t43 = t6 * t31;
t37 = t41 * pkin(3);
t23 = -t37 - pkin(4);
t42 = pkin(4) - t23;
t28 = sin(qJ(5));
t11 = t15 * t28;
t40 = t15 * t31;
t39 = t28 * t31;
t10 = t31 * t16;
t20 = t30 * pkin(3) + qJ(2);
t38 = t15 * t46;
t35 = pkin(4) * t15 - pkin(8) * t16;
t22 = pkin(8) + t44;
t34 = -t15 * t23 - t16 * t22;
t27 = t31 ^ 2;
t26 = t28 ^ 2;
t19 = 0.2e1 * t39;
t9 = t28 * t16;
t8 = t15 * t39;
t7 = t41 * t18 + t29 * t36;
t5 = t6 * t28;
t4 = (t26 - t27) * t15;
t3 = t16 * pkin(4) + t15 * pkin(8) + t20;
t2 = t28 * t3 + t31 * t7;
t1 = -t28 * t7 + t31 * t3;
t12 = [1, 0, 0, -2 * pkin(1), t45, pkin(1) ^ 2 + qJ(2) ^ 2, t32 ^ 2, -0.2e1 * t32 * t30, 0, 0, 0, t30 * t45, t32 * t45, t14, t38, 0, 0, 0, t20 * t46, t20 * t47, t27 * t14, -0.2e1 * t14 * t39, t10 * t47, t28 * t38, t13, 0.2e1 * t1 * t16 - 0.2e1 * t6 * t11, -0.2e1 * t2 * t16 - 0.2e1 * t6 * t40; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48 * t28, t48 * t31; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, t24, -t30 * t33, 0, 0, -t15, -t16, 0, -t6, -t7, -t8, t4, t9, t10, 0, t34 * t28 - t43, t34 * t31 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t40, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t44, t26, t19, 0, 0, 0, -0.2e1 * t23 * t31, 0.2e1 * t23 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, -t6, -t7, -t8, t4, t9, t10, 0, t35 * t28 - t43, t35 * t31 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16, 0, 0, 0, 0, 0, -t40, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t44, t26, t19, 0, 0, 0, t42 * t31, -t42 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t26, t19, 0, 0, 0, 0.2e1 * pkin(4) * t31, -0.2e1 * pkin(4) * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, t11, t16, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * t22, -t31 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t31, 0, -t28 * pkin(8), -t31 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t12;
