% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRP11
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
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:14
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRP11_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP11_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP11_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t30 = sin(qJ(2));
t53 = -0.2e1 * t30;
t32 = cos(qJ(2));
t52 = 0.2e1 * t32;
t51 = 2 * qJ(3);
t33 = -pkin(2) - pkin(7);
t21 = t30 * pkin(6);
t15 = t30 * pkin(3) + t21;
t29 = sin(qJ(4));
t31 = cos(qJ(4));
t37 = -t30 * qJ(3) - pkin(1);
t8 = t33 * t32 + t37;
t5 = t29 * t15 + t31 * t8;
t50 = t30 * pkin(4);
t49 = t29 * t30;
t48 = t29 * t32;
t47 = t29 * t33;
t46 = t30 * t32;
t45 = t30 * t33;
t44 = t31 * t29;
t43 = t31 * t32;
t20 = t31 * t33;
t22 = t32 * pkin(6);
t16 = t32 * pkin(3) + t22;
t24 = t29 ^ 2;
t26 = t31 ^ 2;
t18 = t24 + t26;
t25 = t30 ^ 2;
t27 = t32 ^ 2;
t42 = t25 + t27;
t41 = t30 * qJ(5);
t40 = t32 * qJ(3);
t39 = -0.2e1 * t46;
t38 = -t31 * t15 + t29 * t8;
t2 = t41 + t5;
t3 = t38 - t50;
t1 = t2 * t29 - t3 * t31;
t36 = -t30 * pkin(2) + t40;
t14 = pkin(4) * t31 + t29 * qJ(5);
t35 = t29 * pkin(4) - t31 * qJ(5);
t19 = t31 * t30;
t17 = t30 * t20;
t13 = -t32 * pkin(2) + t37;
t12 = qJ(3) + t35;
t11 = t18 * t33;
t6 = t14 * t32 + t16;
t4 = [1, 0, 0, t25, 0.2e1 * t46, 0, 0, 0, pkin(1) * t52, pkin(1) * t53, 0.2e1 * t42 * pkin(6), t13 * t52, t13 * t53, t42 * pkin(6) ^ 2 + t13 ^ 2, t24 * t27, 0.2e1 * t27 * t44, t29 * t39, t31 * t39, t25, 0.2e1 * t16 * t43 - 0.2e1 * t30 * t38, -0.2e1 * t16 * t48 - 0.2e1 * t5 * t30, -0.2e1 * t3 * t30 + 0.2e1 * t6 * t43, (-t2 * t31 - t29 * t3) * t52, 0.2e1 * t2 * t30 + 0.2e1 * t6 * t48, t2 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, t30, t32, 0, -t21, -t22, t36, t21, t22, t36 * pkin(6), -t29 * t43, (t24 - t26) * t32, t19, -t49, 0, t16 * t29 + t31 * t40 + t17, t16 * t31 + (-t40 - t45) * t29, t12 * t43 + t6 * t29 + t17, -t1, -t6 * t31 + (t12 * t32 + t45) * t29, t1 * t33 + t6 * t12; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(2), t51, pkin(2) ^ 2 + (qJ(3) ^ 2), t26, -0.2e1 * t44, 0, 0, 0, t29 * t51, t31 * t51, 0.2e1 * t12 * t29, -0.2e1 * t11, -0.2e1 * t12 * t31, t18 * t33 ^ 2 + t12 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, 0, 0, t21, 0, 0, 0, 0, 0, t19, -t49, t19, 0, t49, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, -t18, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t43, t30, -t38, -t5, -t38 + 0.2e1 * t50, t35 * t32, 0.2e1 * t41 + t5, -t3 * pkin(4) + t2 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, t20, -t47, t20, -t14, t47, t14 * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, t31, 0, t29, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, -t48, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t4;
