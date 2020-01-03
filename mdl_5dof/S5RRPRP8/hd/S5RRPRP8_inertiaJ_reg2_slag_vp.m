% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP8_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP8_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP8_inertiaJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
t41 = sin(qJ(2));
t38 = t41 ^ 2;
t43 = cos(qJ(2));
t39 = t43 ^ 2;
t54 = t38 + t39;
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t16 = t41 * t40 + t43 * t42;
t53 = 0.2e1 * t16;
t18 = -t43 * t40 + t41 * t42;
t52 = 0.2e1 * t18;
t51 = -0.2e1 * t41;
t50 = -pkin(2) - pkin(3);
t49 = t41 * pkin(6);
t48 = t41 * t43;
t47 = t54 * pkin(6) ^ 2;
t27 = -t43 * pkin(2) - t41 * qJ(3) - pkin(1);
t46 = (pkin(6) - pkin(7)) * t41;
t35 = t43 * pkin(6);
t28 = -t43 * pkin(7) + t35;
t7 = t40 * t28 - t42 * t46;
t11 = t43 * pkin(3) - t27;
t24 = t40 * qJ(3) - t42 * t50;
t45 = -t41 * pkin(2) + t43 * qJ(3);
t9 = t42 * t28 + t40 * t46;
t29 = t40 ^ 2 + t42 ^ 2;
t26 = t42 * qJ(3) + t40 * t50;
t23 = t26 ^ 2;
t22 = 0.2e1 * t54 * pkin(6);
t20 = pkin(4) + t24;
t15 = t18 ^ 2;
t14 = t16 ^ 2;
t13 = 0.2e1 * t26;
t12 = t26 * t40;
t10 = t26 * t16;
t6 = -0.2e1 * t18 * t16;
t5 = t16 * pkin(4) + t11;
t4 = -t40 * t16 - t42 * t18;
t3 = -t16 * qJ(5) + t9;
t1 = t18 * qJ(5) + t7;
t2 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t38, 0.2e1 * t48, 0, t39, 0, 0, 0.2e1 * pkin(1) * t43, pkin(1) * t51, t22, pkin(1) ^ 2 + t47, t38, 0, -0.2e1 * t48, 0, 0, t39, -0.2e1 * t27 * t43, t22, t27 * t51, t27 ^ 2 + t47, t15, t6, 0, t14, 0, 0, t11 * t53, t11 * t52, -0.2e1 * t9 * t16 + 0.2e1 * t7 * t18, t11 ^ 2 + t7 ^ 2 + t9 ^ 2, t15, t6, 0, t14, 0, 0, t5 * t53, t5 * t52, 0.2e1 * t1 * t18 - 0.2e1 * t3 * t16, t1 ^ 2 + t3 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t43, 0, -t49, -t35, 0, 0, 0, t41, 0, 0, -t43, 0, -t49, t45, t35, t45 * pkin(6), 0, 0, -t18, 0, t16, 0, t7, t9, t24 * t18 - t10, t7 * t24 + t9 * t26, 0, 0, -t18, 0, t16, 0, t1, t3, t20 * t18 - t10, t1 * t20 + t3 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2), 0, 0.2e1 * qJ(3), pkin(2) ^ 2 + qJ(3) ^ 2, 0, 0, 0, 0, 0, 1, 0.2e1 * t24, t13, 0, t24 ^ 2 + t23, 0, 0, 0, 0, 0, 1, 0.2e1 * t20, t13, 0, t20 ^ 2 + t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, t49, 0, 0, 0, 0, 0, 0, 0, 0, t4, t9 * t40 - t7 * t42, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t1 * t42 + t3 * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(2), 0, 0, 0, 0, 0, 0, -t42, t40, 0, -t24 * t42 + t12, 0, 0, 0, 0, 0, 0, -t42, t40, 0, -t20 * t42 + t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, 0, -t16, 0, -t7, -t9, 0, 0, 0, 0, t18, 0, -t16, 0, -t1, -t3, -t18 * pkin(4), -t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, -t24, -t26, 0, 0, 0, 0, 0, 0, 0, -1, -0.2e1 * pkin(4) - t24, -t26, 0, -t20 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, 0, 0, 0, 0, 0, 0, 0, 0, t42, -t40, 0, t42 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t2;
