% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x28]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t33 = sin(pkin(9));
t34 = cos(pkin(9));
t37 = sin(qJ(3));
t40 = cos(qJ(3));
t21 = t37 * t33 - t40 * t34;
t22 = t40 * t33 + t37 * t34;
t36 = sin(qJ(4));
t39 = cos(qJ(4));
t14 = t39 * t21 + t36 * t22;
t27 = -t34 * pkin(2) - pkin(1);
t16 = t21 * pkin(3) + t27;
t51 = 0.2e1 * t14 * pkin(4) + 0.2e1 * t16;
t50 = 0.2e1 * t16;
t49 = 0.2e1 * t27;
t35 = sin(qJ(5));
t48 = t35 * pkin(4);
t47 = t36 * pkin(3);
t46 = pkin(6) + qJ(2);
t45 = t33 ^ 2 + t34 ^ 2;
t38 = cos(qJ(5));
t44 = t38 * t47;
t23 = t46 * t33;
t24 = t46 * t34;
t43 = -t40 * t23 - t37 * t24;
t11 = -t22 * pkin(7) + t43;
t42 = t37 * t23 - t40 * t24;
t12 = -t21 * pkin(7) - t42;
t5 = t39 * t11 - t36 * t12;
t30 = t39 * pkin(3);
t28 = t30 + pkin(4);
t17 = t38 * t28 - t35 * t47;
t6 = -t36 * t11 - t39 * t12;
t29 = t38 * pkin(4);
t18 = -t35 * t28 - t44;
t15 = -t36 * t21 + t39 * t22;
t8 = -t35 * t14 + t38 * t15;
t7 = t38 * t14 + t35 * t15;
t4 = -t14 * pkin(8) - t6;
t3 = -t15 * pkin(8) + t5;
t2 = -t35 * t3 - t38 * t4;
t1 = t38 * t3 - t35 * t4;
t9 = [1, 0, 0, 0.2e1 * pkin(1) * t34, -0.2e1 * pkin(1) * t33, 0.2e1 * t45 * qJ(2), t45 * qJ(2) ^ 2 + pkin(1) ^ 2, t22 ^ 2, -0.2e1 * t22 * t21, 0, 0, 0, t21 * t49, t22 * t49, t15 ^ 2, -0.2e1 * t15 * t14, 0, 0, 0, t14 * t50, t15 * t50, t8 ^ 2, -0.2e1 * t8 * t7, 0, 0, 0, t7 * t51, t8 * t51; 0, 0, 0, -t34, t33, 0, -pkin(1), 0, 0, 0, 0, 0, t21, t22, 0, 0, 0, 0, 0, t14, t15, 0, 0, 0, 0, 0, t7, t8; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t22, -t21, 0, t43, t42, 0, 0, t15, -t14, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t30, -0.2e1 * t47, 0, 0, 0, 0, 1, 0.2e1 * t17, 0.2e1 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, t5, t6, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t30, -t47, 0, 0, 0, 0, 1, t17 + t29, -t44 + (-pkin(4) - t28) * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t29, -0.2e1 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, -t7, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t29, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
