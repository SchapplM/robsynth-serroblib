% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPPRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPPRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t19 = (pkin(1) + qJ(3));
t48 = 2 * t19;
t23 = cos(qJ(5));
t47 = 0.2e1 * t23;
t24 = cos(qJ(4));
t46 = 0.2e1 * t24;
t21 = sin(qJ(5));
t45 = t21 * pkin(5);
t44 = t23 * pkin(5);
t14 = t21 ^ 2;
t43 = t14 * t24;
t18 = -pkin(7) + qJ(2);
t42 = t18 * t21;
t22 = sin(qJ(4));
t10 = t21 * t22;
t41 = t21 * t23;
t40 = t21 * t24;
t39 = t22 * t18;
t38 = t23 * t22;
t12 = t23 * t24;
t37 = t24 * t18;
t36 = t24 * t22;
t35 = -qJ(6) - pkin(8);
t16 = t23 ^ 2;
t34 = t14 + t16;
t15 = t22 ^ 2;
t17 = t24 ^ 2;
t33 = -t15 - t17;
t32 = qJ(6) * t24;
t31 = -0.2e1 * t36;
t30 = t18 * t38;
t29 = -pkin(4) * t24 - pkin(8) * t22;
t7 = t22 * pkin(4) - t24 * pkin(8) + t19;
t5 = t23 * t7;
t1 = -t23 * t32 + t5 + (pkin(5) - t42) * t22;
t2 = t30 + (t7 - t32) * t21;
t28 = -t1 * t23 - t2 * t21;
t8 = t35 * t21;
t9 = t35 * t23;
t27 = -t8 * t21 - t9 * t23;
t26 = (qJ(2) ^ 2);
t25 = 2 * qJ(2);
t13 = -pkin(4) - t44;
t11 = t16 * t24;
t6 = (-t18 + t45) * t24;
t4 = t21 * t7 + t30;
t3 = -t21 * t39 + t5;
t20 = [1, 0, 0, -2 * pkin(1), t25, pkin(1) ^ 2 + t26, t25, t48, t19 ^ 2 + t26, t17, t31, 0, 0, 0, t22 * t48, t19 * t46, t16 * t17, -0.2e1 * t17 * t41, t36 * t47, t21 * t31, t15, -0.2e1 * t17 * t42 + 0.2e1 * t3 * t22, -0.2e1 * t17 * t18 * t23 - 0.2e1 * t4 * t22, t28 * t46, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, -1, -t19, 0, 0, 0, 0, 0, -t22, -t24, 0, 0, 0, 0, 0, -t38, t10, t11 + t43, t28; 0, 0, 0, 0, 0, 1, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 1, 0, qJ(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33 * t21, t33 * t23, 0, -t6 * t24 + (-t1 * t21 + t2 * t23) * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t15 + t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, t37, -t39, t21 * t12, t11 - t43, t10, t38, 0, t29 * t21 + t23 * t37, -t21 * t37 + t29 * t23 (-t24 * t8 + t2) * t23 + (t24 * t9 - t1) * t21, t1 * t8 + t6 * t13 - t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 * t9 - t23 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, -t22, 0, 0, 0, 0, 0, t12, -t40, t34 * t22, -t24 * t13 + t27 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t14, 0.2e1 * t41, 0, 0, 0, pkin(4) * t47, -0.2e1 * pkin(4) * t21, 0.2e1 * t27, t13 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t40, t22, t3, -t4, -pkin(5) * t12, t1 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, t21, 0, -t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t38, 0, -pkin(5) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t23, 0, -t21 * pkin(8), -t23 * pkin(8), -t45, t8 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t20;
