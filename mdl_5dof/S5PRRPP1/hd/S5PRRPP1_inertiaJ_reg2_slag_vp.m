% Calculate inertial parameters regressor of joint inertia matrix for
% S5PRRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t24 = sin(pkin(8));
t25 = cos(pkin(8));
t26 = sin(qJ(3));
t27 = cos(qJ(3));
t14 = t24 * t27 + t25 * t26;
t11 = t14 ^ 2;
t12 = t24 * t26 - t25 * t27;
t41 = t12 ^ 2;
t42 = t11 + t41;
t21 = -t27 * pkin(3) - pkin(2);
t40 = 0.2e1 * t21;
t39 = 0.2e1 * t27;
t38 = t24 * pkin(3);
t37 = t25 * pkin(3);
t5 = t12 * t14;
t36 = -qJ(4) - pkin(6);
t22 = t26 ^ 2;
t23 = t27 ^ 2;
t35 = t22 + t23;
t16 = t36 * t27;
t31 = t36 * t26;
t6 = -t24 * t16 - t25 * t31;
t8 = -t25 * t16 + t24 * t31;
t34 = t6 ^ 2 + t8 ^ 2;
t33 = t12 * t6 + t14 * t8;
t30 = -0.2e1 * t8 * t12 + 0.2e1 * t6 * t14;
t19 = pkin(4) + t37;
t17 = qJ(5) + t38;
t3 = t12 * pkin(4) - t14 * qJ(5) + t21;
t1 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t35, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t22, t26 * t39, 0, t23, 0, 0, pkin(2) * t39, -0.2e1 * pkin(2) * t26, 0.2e1 * t35 * pkin(6), t35 * pkin(6) ^ 2 + pkin(2) ^ 2, t11, -0.2e1 * t5, 0, t41, 0, 0, t12 * t40, t14 * t40, t30, t21 ^ 2 + t34, t11, 0, 0.2e1 * t5, 0, 0, t41, 0.2e1 * t3 * t12, t30, -0.2e1 * t3 * t14, t3 ^ 2 + t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, -t26, 0, 0, 0, 0, 0, 0, 0, 0, -t12, -t14, 0, (-t12 * t25 + t14 * t24) * pkin(3), 0, 0, 0, 0, 0, 0, -t12, 0, t14, -t12 * t19 + t14 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, t27, 0, -t26 * pkin(6), -t27 * pkin(6), 0, 0, 0, 0, t14, 0, -t12, 0, -t6, -t8, (-t12 * t24 - t14 * t25) * pkin(3), (t24 * t8 - t25 * t6) * pkin(3), 0, t14, 0, 0, t12, 0, -t6, -t17 * t12 - t19 * t14, t8, t8 * t17 - t6 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t38, 0, (t24 ^ 2 + t25 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t19, 0, 0.2e1 * t17, t17 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, t14, 0, t21, 0, 0, 0, 0, 0, 0, t12, 0, -t14, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
