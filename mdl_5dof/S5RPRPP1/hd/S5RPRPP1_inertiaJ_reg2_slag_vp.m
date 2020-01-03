% Calculate inertial parameters regressor of joint inertia matrix for
% S5RPRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRPP1_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t27 = sin(pkin(8));
t29 = cos(pkin(8));
t31 = sin(qJ(3));
t32 = cos(qJ(3));
t16 = t27 * t32 + t29 * t31;
t13 = t16 ^ 2;
t14 = t27 * t31 - t29 * t32;
t48 = t14 ^ 2;
t49 = t13 + t48;
t30 = cos(pkin(7));
t42 = t30 * pkin(1);
t24 = -pkin(2) - t42;
t17 = -t32 * pkin(3) + t24;
t47 = 0.2e1 * t17;
t46 = 0.2e1 * t31;
t45 = t27 * pkin(3);
t28 = sin(pkin(7));
t44 = t28 * pkin(1);
t43 = t29 * pkin(3);
t8 = t16 * t14;
t25 = t31 ^ 2;
t26 = t32 ^ 2;
t41 = t25 + t26;
t21 = pkin(6) + t44;
t40 = qJ(4) + t21;
t10 = t40 * t32;
t36 = t40 * t31;
t5 = t27 * t10 + t29 * t36;
t7 = t29 * t10 - t27 * t36;
t39 = t5 ^ 2 + t7 ^ 2;
t38 = t5 * t14 + t7 * t16;
t35 = -0.2e1 * t7 * t14 + 0.2e1 * t5 * t16;
t22 = pkin(4) + t43;
t19 = qJ(5) + t45;
t3 = t14 * pkin(4) - t16 * qJ(5) + t17;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t42, -0.2e1 * t44, 0, (t28 ^ 2 + t30 ^ 2) * pkin(1) ^ 2, t25, t32 * t46, 0, t26, 0, 0, -0.2e1 * t24 * t32, t24 * t46, 0.2e1 * t41 * t21, t41 * t21 ^ 2 + t24 ^ 2, t13, -0.2e1 * t8, 0, t48, 0, 0, t14 * t47, t16 * t47, t35, t17 ^ 2 + t39, t13, 0, 0.2e1 * t8, 0, 0, t48, 0.2e1 * t3 * t14, t35, -0.2e1 * t3 * t16, t3 ^ 2 + t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, t32, 0, -t31 * t21, -t32 * t21, 0, 0, 0, 0, t16, 0, -t14, 0, -t5, -t7, (-t14 * t27 - t16 * t29) * pkin(3), (t27 * t7 - t29 * t5) * pkin(3), 0, t16, 0, 0, t14, 0, -t5, -t19 * t14 - t22 * t16, t7, t7 * t19 - t5 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, -t31, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t16, 0, (-t14 * t29 + t16 * t27) * pkin(3), 0, 0, 0, 0, 0, 0, -t14, 0, t16, -t14 * t22 + t16 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t43, -0.2e1 * t45, 0, (t27 ^ 2 + t29 ^ 2) * pkin(3) ^ 2, 0, 0, 0, 1, 0, 0, 0.2e1 * t22, 0, 0.2e1 * t19, t19 ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t16, 0, t17, 0, 0, 0, 0, 0, 0, t14, 0, -t16, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t1;
