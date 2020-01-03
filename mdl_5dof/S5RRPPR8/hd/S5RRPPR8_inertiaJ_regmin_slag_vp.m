% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPPR8_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR8_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR8_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t39 = sin(qJ(2));
t41 = cos(qJ(2));
t26 = -t41 * pkin(2) - t39 * qJ(3) - pkin(1);
t14 = t41 * pkin(3) - t26;
t36 = sin(pkin(8));
t37 = cos(pkin(8));
t16 = t39 * t36 + t41 * t37;
t51 = 0.2e1 * t16 * pkin(4) + 0.2e1 * t14;
t50 = 0.2e1 * t14;
t49 = -0.2e1 * t39;
t48 = 0.2e1 * t41;
t47 = -pkin(2) - pkin(3);
t46 = t39 * pkin(6);
t31 = t41 * pkin(6);
t27 = -t41 * qJ(4) + t31;
t44 = (pkin(6) - qJ(4)) * t39;
t12 = t37 * t27 + t36 * t44;
t34 = t39 ^ 2;
t45 = t41 ^ 2 + t34;
t10 = t36 * t27 - t37 * t44;
t23 = t36 * qJ(3) - t37 * t47;
t43 = -t39 * pkin(2) + t41 * qJ(3);
t40 = cos(qJ(5));
t38 = sin(qJ(5));
t25 = t37 * qJ(3) + t36 * t47;
t22 = -pkin(4) - t23;
t18 = t40 * t36 + t38 * t37;
t17 = -t41 * t36 + t39 * t37;
t15 = t38 * t36 - t40 * t37;
t8 = t38 * t22 + t40 * t25;
t7 = -t40 * t22 + t38 * t25;
t6 = -t38 * t16 + t40 * t17;
t5 = t40 * t16 + t38 * t17;
t4 = -t16 * pkin(7) + t12;
t3 = -t17 * pkin(7) - t10;
t2 = t38 * t3 + t40 * t4;
t1 = -t40 * t3 + t38 * t4;
t9 = [1, 0, 0, t34, t39 * t48, 0, 0, 0, pkin(1) * t48, pkin(1) * t49, -0.2e1 * t26 * t41, 0.2e1 * t45 * pkin(6), t26 * t49, t45 * pkin(6) ^ 2 + t26 ^ 2, t16 * t50, t17 * t50, 0.2e1 * t10 * t17 - 0.2e1 * t12 * t16, t10 ^ 2 + t12 ^ 2 + t14 ^ 2, t6 ^ 2, -0.2e1 * t6 * t5, 0, 0, 0, t5 * t51, t6 * t51; 0, 0, 0, 0, 0, t39, t41, 0, -t46, -t31, -t46, t43, t31, t43 * pkin(6), t10, t12, -t25 * t16 + t23 * t17, t10 * t23 + t12 * t25, 0, 0, -t6, t5, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0.2e1 * t23, 0.2e1 * t25, 0, t23 ^ 2 + t25 ^ 2, 0, 0, 0, 0, 1, 0.2e1 * t7, 0.2e1 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, 0, t46, 0, 0, -t36 * t16 - t37 * t17, -t10 * t37 + t12 * t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), -t37, t36, 0, -t23 * t37 + t25 * t36, 0, 0, 0, 0, 0, t15, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t36 ^ 2 + t37 ^ 2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t17, 0, t14, 0, 0, 0, 0, 0, t5, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t7, -t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t15, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
