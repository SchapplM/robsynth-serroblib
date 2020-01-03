% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t36 = cos(qJ(4));
t43 = t36 * pkin(4);
t28 = cos(qJ(2)) * pkin(1);
t26 = t28 + pkin(2);
t30 = sin(pkin(9));
t31 = cos(pkin(9));
t45 = sin(qJ(2)) * pkin(1);
t11 = t31 * t26 - t30 * t45;
t9 = -pkin(3) - t11;
t8 = t9 - t43;
t50 = 0.2e1 * t8;
t24 = -t31 * pkin(2) - pkin(3);
t18 = t24 - t43;
t49 = 0.2e1 * t18;
t33 = sin(qJ(4));
t48 = 0.2e1 * t33;
t47 = -0.2e1 * t36;
t32 = sin(qJ(5));
t46 = t32 * pkin(4);
t35 = cos(qJ(5));
t44 = t35 * pkin(4);
t42 = t18 + t8;
t41 = t24 + t9;
t12 = t30 * t26 + t31 * t45;
t10 = pkin(7) + t12;
t40 = t36 * t10;
t23 = t30 * pkin(2) + pkin(7);
t39 = t36 * t23;
t29 = t33 ^ 2;
t27 = t36 * pkin(8);
t22 = t36 * t48;
t17 = t32 * t36 + t35 * t33;
t16 = t32 * t33 - t35 * t36;
t15 = t17 ^ 2;
t14 = t27 + t39;
t13 = (-pkin(8) - t23) * t33;
t7 = t27 + t40;
t6 = (-pkin(8) - t10) * t33;
t5 = -0.2e1 * t17 * t16;
t4 = -t32 * t13 - t35 * t14;
t3 = t35 * t13 - t32 * t14;
t2 = -t32 * t6 - t35 * t7;
t1 = -t32 * t7 + t35 * t6;
t19 = [1, 0, 0, 1, 0.2e1 * t28, -0.2e1 * t45, t11 ^ 2 + t12 ^ 2, t29, t22, 0, 0, 0, t9 * t47, t9 * t48, t15, t5, 0, 0, 0, t16 * t50, t17 * t50; 0, 0, 0, 1, t28, -t45, (t11 * t31 + t12 * t30) * pkin(2), t29, t22, 0, 0, 0, -t41 * t36, t41 * t33, t15, t5, 0, 0, 0, t42 * t16, t42 * t17; 0, 0, 0, 1, 0, 0, (t30 ^ 2 + t31 ^ 2) * pkin(2) ^ 2, t29, t22, 0, 0, 0, t24 * t47, t24 * t48, t15, t5, 0, 0, 0, t16 * t49, t17 * t49; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36, 0, -t33 * t10, -t40, 0, 0, t17, -t16, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t36, 0, -t33 * t23, -t39, 0, 0, t17, -t16, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t33, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t44, -0.2e1 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t16, -t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t44, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t19;
