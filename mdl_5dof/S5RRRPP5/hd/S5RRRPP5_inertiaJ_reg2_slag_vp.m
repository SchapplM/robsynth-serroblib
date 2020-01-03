% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRPP5_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t33 = sin(qJ(3));
t34 = sin(qJ(2));
t35 = cos(qJ(3));
t36 = cos(qJ(2));
t14 = t33 * t34 - t35 * t36;
t12 = t14 ^ 2;
t29 = -t36 * pkin(2) - pkin(1);
t56 = 0.2e1 * t29;
t55 = 0.2e1 * t36;
t37 = pkin(3) + pkin(4);
t54 = -pkin(7) - pkin(6);
t53 = t35 * pkin(2);
t16 = t33 * t36 + t35 * t34;
t52 = t16 * t14;
t30 = t33 * pkin(2);
t24 = t30 + qJ(4);
t51 = t24 * t14;
t31 = t34 ^ 2;
t32 = t36 ^ 2;
t50 = t31 + t32;
t49 = qJ(4) * t14;
t18 = t54 * t36;
t47 = t54 * t34;
t10 = -t35 * t18 + t33 * t47;
t8 = -t33 * t18 - t35 * t47;
t48 = t10 ^ 2 + t8 ^ 2;
t27 = pkin(3) + t53;
t41 = 0.2e1 * pkin(3);
t46 = t41 + t53;
t45 = -0.2e1 * t10 * t14 + 0.2e1 * t8 * t16;
t44 = t16 * qJ(4) - t29;
t40 = qJ(4) ^ 2;
t39 = 0.2e1 * qJ(4);
t25 = t39 + t30;
t22 = pkin(4) + t27;
t21 = t24 ^ 2;
t20 = t24 * qJ(4);
t19 = 0.2e1 * t24;
t13 = t16 ^ 2;
t6 = 0.2e1 * t52;
t5 = t14 * pkin(3) - t44;
t3 = t14 * qJ(5) + t10;
t2 = -t16 * qJ(5) + t8;
t1 = -t37 * t14 + t44;
t4 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t31, t34 * t55, 0, t32, 0, 0, pkin(1) * t55, -0.2e1 * pkin(1) * t34, 0.2e1 * t50 * pkin(6), t50 * pkin(6) ^ 2 + pkin(1) ^ 2, t13, -0.2e1 * t52, 0, t12, 0, 0, t14 * t56, t16 * t56, t45, t29 ^ 2 + t48, t13, 0, t6, 0, 0, t12, 0.2e1 * t5 * t14, t45, -0.2e1 * t5 * t16, t5 ^ 2 + t48, t13, t6, 0, t12, 0, 0, -0.2e1 * t1 * t14, 0.2e1 * t1 * t16, 0.2e1 * t3 * t14 - 0.2e1 * t2 * t16, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, t36, 0, -t34 * pkin(6), -t36 * pkin(6), 0, 0, 0, 0, t16, 0, -t14, 0, -t8, -t10, (-t14 * t33 - t16 * t35) * pkin(2), (t10 * t33 - t35 * t8) * pkin(2), 0, t16, 0, 0, t14, 0, -t8, -t27 * t16 - t51, t10, t10 * t24 - t8 * t27, 0, 0, -t16, 0, -t14, 0, -t2, t3, t22 * t16 + t51, -t2 * t22 + t3 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t53, -0.2e1 * t30, 0, (t33 ^ 2 + t35 ^ 2) * pkin(2) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t27, 0, t19, t27 ^ 2 + t21, 0, 0, 0, 0, 0, 1, 0.2e1 * t22, t19, 0, t22 ^ 2 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, -t14, 0, -t8, -t10, 0, 0, 0, t16, 0, 0, t14, 0, -t8, -pkin(3) * t16 - t49, t10, -t8 * pkin(3) + t10 * qJ(4), 0, 0, -t16, 0, -t14, 0, -t2, t3, t37 * t16 + t49, t3 * qJ(4) - t2 * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t53, -t30, 0, 0, 0, 0, 0, 1, 0, 0, t46, 0, t25, t27 * pkin(3) + t20, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4) + t46, t25, 0, t22 * t37 + t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t41, 0, t39, pkin(3) ^ 2 + t40, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, t39, 0, t37 ^ 2 + t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t16, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t27, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(3), 0, 0, 0, 0, 0, 0, -1, 0, 0, -t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, t16, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
