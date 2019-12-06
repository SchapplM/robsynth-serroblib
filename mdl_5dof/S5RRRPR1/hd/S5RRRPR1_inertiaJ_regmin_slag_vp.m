% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x26]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:39
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t44 = sin(qJ(3));
t45 = sin(qJ(2));
t47 = cos(qJ(2));
t50 = cos(qJ(3));
t29 = t44 * t45 - t50 * t47;
t30 = t44 * t47 + t50 * t45;
t41 = sin(pkin(9));
t42 = cos(pkin(9));
t18 = -t42 * t29 - t41 * t30;
t39 = -t47 * pkin(2) - pkin(1);
t23 = t29 * pkin(3) + t39;
t56 = -0.2e1 * t18 * pkin(4) + 0.2e1 * t23;
t55 = 0.2e1 * t39;
t54 = 0.2e1 * t47;
t53 = pkin(6) + pkin(7);
t52 = pkin(3) * t41;
t51 = t44 * pkin(2);
t34 = t53 * t45;
t35 = t53 * t47;
t20 = -t50 * t34 - t44 * t35;
t15 = -t30 * qJ(4) + t20;
t21 = t44 * t34 - t50 * t35;
t16 = -t29 * qJ(4) - t21;
t6 = t41 * t15 + t42 * t16;
t40 = t50 * pkin(2);
t38 = t40 + pkin(3);
t27 = t41 * t38 + t42 * t51;
t49 = -t27 - t52;
t5 = t42 * t15 - t41 * t16;
t25 = t42 * t38 - t41 * t51;
t46 = cos(qJ(5));
t43 = sin(qJ(5));
t36 = t42 * pkin(3) + pkin(4);
t33 = t46 * t36;
t28 = -t43 * t36 - t46 * t52;
t26 = -t43 * t52 + t33;
t24 = pkin(4) + t25;
t22 = t46 * t24;
t19 = -t41 * t29 + t42 * t30;
t14 = -t43 * t24 - t46 * t27;
t13 = -t43 * t27 + t22;
t8 = t43 * t18 + t46 * t19;
t7 = -t46 * t18 + t43 * t19;
t4 = t18 * pkin(8) + t6;
t3 = -t19 * pkin(8) + t5;
t2 = -t43 * t3 - t46 * t4;
t1 = t46 * t3 - t43 * t4;
t9 = [1, 0, 0, t45 ^ 2, t45 * t54, 0, 0, 0, pkin(1) * t54, -0.2e1 * pkin(1) * t45, t30 ^ 2, -0.2e1 * t30 * t29, 0, 0, 0, t29 * t55, t30 * t55, 0.2e1 * t6 * t18 - 0.2e1 * t5 * t19, t23 ^ 2 + t5 ^ 2 + t6 ^ 2, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t56, t8 * t56; 0, 0, 0, 0, 0, t45, t47, 0, -t45 * pkin(6), -t47 * pkin(6), 0, 0, t30, -t29, 0, t20, t21, t27 * t18 - t25 * t19, t5 * t25 + t6 * t27, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t40, -0.2e1 * t51, 0, t25 ^ 2 + t27 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t13, 0.2e1 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, 0, t20, t21, (t18 * t41 - t19 * t42) * pkin(3), (t41 * t6 + t42 * t5) * pkin(3), 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t40, -t51, 0, (t25 * t42 + t27 * t41) * pkin(3), 0, 0, 0, 0, 1, t49 * t43 + t22 + t33, t49 * t46 + (-t24 - t36) * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t41 ^ 2 + t42 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t26, 0.2e1 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
