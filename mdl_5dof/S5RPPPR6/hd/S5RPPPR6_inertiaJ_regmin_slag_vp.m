% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPPR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t31 = cos(pkin(7));
t27 = t31 ^ 2;
t29 = sin(pkin(7));
t51 = t29 ^ 2 + t27;
t50 = -0.2e1 * t29;
t49 = 0.2e1 * t31;
t36 = -t29 * qJ(3) - pkin(1);
t11 = (-pkin(2) - qJ(4)) * t31 + t36;
t20 = t29 * qJ(2);
t15 = t29 * pkin(3) + t20;
t28 = sin(pkin(8));
t30 = cos(pkin(8));
t6 = t30 * t11 + t28 * t15;
t26 = t30 ^ 2;
t48 = t26 * t31;
t47 = t28 * t29;
t46 = t28 * t31;
t45 = t30 * t29;
t18 = t30 * t31;
t32 = sin(qJ(5));
t44 = t32 * t28;
t43 = t32 * t30;
t33 = cos(qJ(5));
t42 = t33 * t28;
t41 = t33 * t30;
t40 = t51 * qJ(2) ^ 2;
t16 = (pkin(3) + qJ(2)) * t31;
t17 = t28 ^ 2 + t26;
t39 = 0.2e1 * t18;
t38 = t31 * t44;
t37 = t31 * t42;
t5 = -t28 * t11 + t30 * t15;
t35 = t5 * t28 - t6 * t30;
t14 = -t31 * pkin(2) + t36;
t13 = 0.2e1 * t51 * qJ(2);
t10 = t29 * t32 - t37;
t9 = t29 * t33 + t38;
t7 = (pkin(4) * t30 + pkin(6) * t28) * t31 + t16;
t4 = t29 * pkin(6) + t6;
t3 = -t29 * pkin(4) - t5;
t2 = t32 * t7 + t33 * t4;
t1 = -t32 * t4 + t33 * t7;
t8 = [1, 0, 0, pkin(1) * t49, pkin(1) * t50, t13, pkin(1) ^ 2 + t40, t13, t14 * t49, t14 * t50, t14 ^ 2 + t40, 0.2e1 * t16 * t18 + 0.2e1 * t5 * t29, -0.2e1 * t16 * t46 - 0.2e1 * t6 * t29, t35 * t49, t16 ^ 2 + t5 ^ 2 + t6 ^ 2, t10 ^ 2, 0.2e1 * t10 * t9, t10 * t39, t9 * t39, t26 * t27, 0.2e1 * t1 * t18 - 0.2e1 * t3 * t9, 0.2e1 * t3 * t10 - 0.2e1 * t2 * t18; 0, 0, 0, -t31, t29, 0, -pkin(1), 0, t31, -t29, t14, -t47, -t45, -t17 * t31, -t35, 0, 0, 0, 0, 0, -t28 * t9 - t32 * t48, t28 * t10 - t33 * t48; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t29, 0, 0, t20, t45, -t47, 0, t6 * t28 + t5 * t30, 0, 0, 0, 0, 0, (t9 - t38) * t30, (-t10 - t37) * t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, t17, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t46, 0, t16, 0, 0, 0, 0, 0, t31 * t41, -t31 * t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, t9, t18, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t41; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t8;
