% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:30
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t29 = sin(pkin(8));
t30 = cos(pkin(8));
t32 = sin(qJ(2));
t34 = cos(qJ(2));
t17 = t29 * t32 - t30 * t34;
t18 = t29 * t34 + t30 * t32;
t28 = -t34 * pkin(2) - pkin(1);
t36 = t18 * qJ(4) - t28;
t42 = 0.2e1 * (-pkin(3) - pkin(4)) * t17 + 0.2e1 * t36;
t41 = 0.2e1 * t34;
t40 = -qJ(3) - pkin(6);
t22 = t40 * t34;
t38 = t40 * t32;
t13 = -t30 * t22 + t29 * t38;
t11 = -t29 * t22 - t30 * t38;
t39 = t11 ^ 2 + t13 ^ 2;
t26 = t30 * pkin(2) + pkin(3);
t37 = 0.2e1 * t11 * t18 - 0.2e1 * t13 * t17;
t33 = cos(qJ(5));
t31 = sin(qJ(5));
t24 = t29 * pkin(2) + qJ(4);
t23 = -pkin(4) - t26;
t15 = t31 * t23 + t33 * t24;
t14 = -t33 * t23 + t31 * t24;
t9 = t31 * t17 + t33 * t18;
t8 = -t33 * t17 + t31 * t18;
t7 = t17 * pkin(3) - t36;
t5 = t17 * pkin(7) + t13;
t4 = -t18 * pkin(7) + t11;
t2 = t31 * t4 + t33 * t5;
t1 = t31 * t5 - t33 * t4;
t3 = [1, 0, 0, t32 ^ 2, t32 * t41, 0, 0, 0, pkin(1) * t41, -0.2e1 * pkin(1) * t32, t37, t28 ^ 2 + t39, 0.2e1 * t7 * t17, t37, -0.2e1 * t7 * t18, t7 ^ 2 + t39, t9 ^ 2, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t42, t9 * t42; 0, 0, 0, 0, 0, t32, t34, 0, -t32 * pkin(6), -t34 * pkin(6), (-t17 * t29 - t18 * t30) * pkin(2), (-t11 * t30 + t13 * t29) * pkin(2), -t11, -t24 * t17 - t26 * t18, t13, -t11 * t26 + t13 * t24, 0, 0, -t9, t8, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, (t29 ^ 2 + t30 ^ 2) * pkin(2) ^ 2, 0.2e1 * t26, 0, 0.2e1 * t24, t24 ^ 2 + t26 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t14, 0.2e1 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, t17, 0, -t18, t7, 0, 0, 0, 0, 0, -t8, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, t11, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t26, 0, 0, 0, 0, 0, -t33, t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t14, -t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t3;
