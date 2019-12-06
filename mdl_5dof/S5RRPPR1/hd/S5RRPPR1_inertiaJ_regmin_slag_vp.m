% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(pkin(8));
t24 = t33 * pkin(2) + qJ(4);
t32 = sin(pkin(9));
t34 = cos(pkin(9));
t41 = t32 ^ 2 + t34 ^ 2;
t42 = t41 * t24;
t46 = t34 * pkin(4);
t29 = cos(qJ(2)) * pkin(1);
t27 = t29 + pkin(2);
t35 = cos(pkin(8));
t45 = sin(qJ(2)) * pkin(1);
t10 = t35 * t27 - t33 * t45;
t9 = -pkin(3) - t10;
t4 = t9 - t46;
t51 = 0.2e1 * t4;
t26 = -t35 * pkin(2) - pkin(3);
t17 = t26 - t46;
t50 = 0.2e1 * t17;
t49 = 0.2e1 * t32;
t48 = -0.2e1 * t34;
t11 = t33 * t27 + t35 * t45;
t8 = qJ(4) + t11;
t47 = t41 * t8;
t44 = t17 + t4;
t43 = t26 + t9;
t38 = cos(qJ(5));
t36 = sin(qJ(5));
t28 = t34 * pkin(7);
t16 = t38 * t32 + t36 * t34;
t15 = t36 * t32 - t38 * t34;
t14 = t16 ^ 2;
t13 = t34 * t24 + t28;
t12 = (-pkin(7) - t24) * t32;
t3 = t34 * t8 + t28;
t2 = (-pkin(7) - t8) * t32;
t1 = -0.2e1 * t16 * t15;
t5 = [1, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t45, t10 ^ 2 + t11 ^ 2, t9 * t48, t9 * t49, 0.2e1 * t47, t41 * t8 ^ 2 + t9 ^ 2, t14, t1, 0, 0, 0, t15 * t51, t16 * t51; 0, 0, 0, 1, t29, -t45, (t10 * t35 + t11 * t33) * pkin(2), -t43 * t34, t43 * t32, t42 + t47, t9 * t26 + t8 * t42, t14, t1, 0, 0, 0, t44 * t15, t44 * t16; 0, 0, 0, 1, 0, 0, (t33 ^ 2 + t35 ^ 2) * pkin(2) ^ 2, t26 * t48, t26 * t49, 0.2e1 * t42, t41 * t24 ^ 2 + t26 ^ 2, t14, t1, 0, 0, 0, t15 * t50, t16 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, -t34, t32, 0, t9, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, -t34, t32, 0, t26, 0, 0, 0, 0, 0, t15, t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t38 * t2 - t36 * t3, -t36 * t2 - t38 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t38 * t12 - t36 * t13, -t36 * t12 - t38 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t5;
