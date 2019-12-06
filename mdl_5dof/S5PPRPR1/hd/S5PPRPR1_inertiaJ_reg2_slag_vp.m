% Calculate inertial parameters regressor of joint inertia matrix for
% S5PPRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRPR1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRPR1_inertiaJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t23 = sin(pkin(8));
t25 = cos(pkin(8));
t27 = sin(qJ(3));
t34 = cos(qJ(3));
t11 = t27 * t23 - t34 * t25;
t6 = t11 ^ 2;
t24 = cos(pkin(9));
t19 = -t24 * pkin(4) - pkin(3);
t36 = 0.2e1 * t19;
t35 = 0.2e1 * t24;
t33 = cos(qJ(5));
t32 = pkin(6) + qJ(4);
t22 = sin(pkin(9));
t20 = t22 ^ 2;
t21 = t24 ^ 2;
t31 = t20 + t21;
t30 = t31 * qJ(4);
t26 = sin(qJ(5));
t13 = t33 * t22 + t26 * t24;
t29 = -t26 * t22 + t33 * t24;
t16 = t32 * t24;
t15 = t32 * t22;
t14 = t34 * t23 + t27 * t25;
t8 = t14 ^ 2;
t7 = t13 ^ 2;
t5 = t29 ^ 2;
t4 = -t26 * t15 + t33 * t16;
t3 = -t33 * t15 - t26 * t16;
t2 = t29 * t14;
t1 = t13 * t14;
t9 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 + t25 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8 + t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31 * t8 + t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1 * t29 + t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7 + t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t14, 0, 0, 0, 0, 0, 0, 0, 0, -t11 * t24, t11 * t22, t31 * t14, -t11 * pkin(3) + t14 * t30, 0, 0, 0, 0, 0, 0, -t11 * t29, t11 * t13, t1 * t13 + t2 * t29, -t1 * t3 + t11 * t19 + t2 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13 * t4 + t29 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t20, t22 * t35, 0, t21, 0, 0, pkin(3) * t35, -0.2e1 * pkin(3) * t22, 0.2e1 * t30, t31 * qJ(4) ^ 2 + pkin(3) ^ 2, t7, 0.2e1 * t13 * t29, 0, t5, 0, 0, -t29 * t36, t13 * t36, -0.2e1 * t3 * t13 + 0.2e1 * t29 * t4, t19 ^ 2 + t3 ^ 2 + t4 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t24, t22, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -t29, t13, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t13, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, 0, t29, 0, t3, -t4, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t9;
