% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:18
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t28 = sin(pkin(9));
t48 = -0.2e1 * t28;
t47 = 0.2e1 * t28;
t45 = sin(qJ(3)) * pkin(2);
t22 = qJ(4) + t45;
t29 = cos(pkin(9));
t30 = sin(qJ(5));
t32 = cos(qJ(5));
t39 = t32 * t29;
t10 = -t29 * pkin(4) - t28 * pkin(7) - pkin(3);
t44 = cos(qJ(3)) * pkin(2);
t7 = t10 - t44;
t3 = t22 * t39 + t30 * t7;
t26 = t28 ^ 2;
t41 = t26 * t32;
t46 = t22 * t41 + t3 * t29;
t23 = -pkin(3) - t44;
t43 = pkin(3) - t23;
t24 = t26 * qJ(4);
t35 = qJ(4) * t29;
t6 = t30 * t10 + t32 * t35;
t42 = t32 * t24 + t6 * t29;
t13 = t26 * t22;
t40 = t30 * t28;
t19 = t30 * t29;
t27 = t29 ^ 2;
t38 = t27 * t22 + t13;
t37 = t27 * qJ(4) + t24;
t36 = t26 + t27;
t21 = t32 * t28;
t20 = t32 ^ 2 * t26;
t16 = t30 * t24;
t15 = -0.2e1 * t30 * t41;
t12 = t39 * t48;
t11 = t19 * t47;
t8 = t30 * t13;
t5 = t32 * t10 - t30 * t35;
t2 = -t22 * t19 + t32 * t7;
t1 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 1, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t45, -0.2e1 * t23 * t29, t23 * t47, 0.2e1 * t38, t36 * t22 ^ 2 + t23 ^ 2, t20, t15, t12, t11, t27, -0.2e1 * t2 * t29 + 0.2e1 * t8, 0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 1, t44, -t45, t43 * t29, -t43 * t28, t37 + t38, t36 * t22 * qJ(4) - t23 * pkin(3), t20, t15, t12, t11, t27, t16 + t8 + (-t2 - t5) * t29, t42 + t46; 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(3) * t29, pkin(3) * t48, 0.2e1 * t37, t36 * qJ(4) ^ 2 + pkin(3) ^ 2, t20, t15, t12, t11, t27, -0.2e1 * t5 * t29 + 0.2e1 * t16, 0.2e1 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, t23, 0, 0, 0, 0, 0, -t39, t19; 0, 0, 0, 0, 0, 0, 0, -t29, t28, 0, -pkin(3), 0, 0, 0, 0, 0, -t39, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t40, -t29, t2, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t40, -t29, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t1;
