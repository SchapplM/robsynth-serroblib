% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x18]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t21 = sin(qJ(3));
t38 = t21 * pkin(2);
t14 = pkin(7) + t38;
t20 = sin(qJ(4));
t18 = t20 ^ 2;
t23 = cos(qJ(4));
t28 = t23 ^ 2 + t18;
t45 = t28 * t14;
t44 = -0.2e1 * t20;
t43 = -0.2e1 * t23;
t42 = 0.2e1 * t23;
t24 = cos(qJ(3));
t35 = t24 * pkin(2);
t7 = -t23 * pkin(4) - t20 * qJ(5) - pkin(3);
t3 = t7 - t35;
t41 = -t3 - t7;
t40 = t20 * pkin(7);
t22 = sin(qJ(2));
t25 = cos(qJ(2));
t6 = t21 * t25 + t24 * t22;
t39 = t20 * t6;
t37 = t23 * pkin(7);
t36 = t23 * t6;
t5 = t21 * t22 - t24 * t25;
t34 = t5 * t23;
t15 = -pkin(3) - t35;
t33 = pkin(3) - t15;
t31 = t20 * t14;
t30 = t23 * t14;
t29 = t28 * pkin(7);
t1 = t28 * t6;
t8 = -t20 * pkin(4) + t23 * qJ(5);
t11 = t20 * t42;
t2 = t5 * t20;
t4 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28 * t6 ^ 2 + t5 ^ 2; 0, 0, t25, -t22, 0, -t5, -t6, 0, 0, 0, 0, 0, -t34, t2, -t34, t1, -t2, t5 * t3 + t45 * t6; 0, 1, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t38, t18, t11, 0, 0, 0, t15 * t43, 0.2e1 * t15 * t20, t3 * t43, 0.2e1 * t45, t3 * t44, t28 * t14 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, -t5, -t6, 0, 0, 0, 0, 0, -t34, t2, -t34, t1, -t2, pkin(7) * t1 + t5 * t7; 0, 0, 0, 0, 1, t35, -t38, t18, t11, 0, 0, 0, t33 * t23, -t33 * t20, t41 * t23, t29 + t45, t41 * t20, pkin(7) * t45 + t3 * t7; 0, 0, 0, 0, 1, 0, 0, t18, t11, 0, 0, 0, pkin(3) * t42, pkin(3) * t44, t7 * t43, 0.2e1 * t29, t7 * t44, t28 * pkin(7) ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t39, -t36, -t39, 0, t36, t8 * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t23, 0, -t31, -t30, -t31, t8, t30, t8 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, t23, 0, -t40, -t37, -t40, t8, t37, t8 * pkin(7); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, 0, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
