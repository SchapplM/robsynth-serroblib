% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x25]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRR6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:01:40
% EndTime: 2019-12-31 19:01:41
% DurationCPUTime: 0.28s
% Computational Cost: add. (200->50), mult. (405->88), div. (0->0), fcn. (481->8), ass. (0->46)
t26 = cos(pkin(9));
t19 = -t26 * pkin(1) - pkin(2);
t31 = cos(qJ(3));
t16 = -t31 * pkin(3) + t19;
t48 = 0.2e1 * t16;
t29 = sin(qJ(3));
t47 = 0.2e1 * t29;
t28 = sin(qJ(4));
t46 = t28 * pkin(3);
t30 = cos(qJ(5));
t25 = sin(pkin(9));
t18 = t25 * pkin(1) + pkin(6);
t43 = pkin(7) + t18;
t12 = t43 * t31;
t42 = cos(qJ(4));
t35 = t42 * t29;
t5 = t28 * t12 + t43 * t35;
t45 = t5 * t30;
t36 = t42 * pkin(3);
t22 = -t36 - pkin(4);
t44 = pkin(4) - t22;
t15 = t28 * t31 + t35;
t27 = sin(qJ(5));
t41 = t27 * t15;
t40 = t27 * t30;
t39 = t28 * t29;
t38 = t30 * t15;
t14 = -t42 * t31 + t39;
t37 = -0.2e1 * t15 * t14;
t34 = -pkin(4) * t15 - pkin(8) * t14;
t21 = pkin(8) + t46;
t33 = -t14 * t21 + t15 * t22;
t24 = t30 ^ 2;
t23 = t27 ^ 2;
t17 = 0.2e1 * t40;
t13 = t15 ^ 2;
t11 = t30 * t14;
t10 = t27 * t14;
t8 = t27 * t38;
t7 = (-t23 + t24) * t15;
t6 = t42 * t12 - t43 * t39;
t4 = t14 * pkin(4) - t15 * pkin(8) + t16;
t3 = t5 * t27;
t2 = t27 * t4 + t30 * t6;
t1 = -t27 * t6 + t30 * t4;
t9 = [1, 0, 0, (t25 ^ 2 + t26 ^ 2) * pkin(1) ^ 2, t29 ^ 2, t31 * t47, 0, 0, 0, -0.2e1 * t19 * t31, t19 * t47, t13, t37, 0, 0, 0, t14 * t48, t15 * t48, t24 * t13, -0.2e1 * t13 * t40, 0.2e1 * t14 * t38, t27 * t37, t14 ^ 2, 0.2e1 * t1 * t14 + 0.2e1 * t5 * t41, -0.2e1 * t2 * t14 + 0.2e1 * t5 * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, t29, t31, 0, -t29 * t18, -t31 * t18, 0, 0, t15, -t14, 0, -t5, -t6, t8, t7, t10, t11, 0, t33 * t27 - t45, t33 * t30 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t29, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t36, -0.2e1 * t46, t23, t17, 0, 0, 0, -0.2e1 * t22 * t30, 0.2e1 * t22 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t15, -t14, 0, -t5, -t6, t8, t7, t10, t11, 0, t34 * t27 - t45, t34 * t30 + t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14, -t15, 0, 0, 0, 0, 0, -t11, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t36, -t46, t23, t17, 0, 0, 0, t44 * t30, -t44 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, t17, 0, 0, 0, 0.2e1 * pkin(4) * t30, -0.2e1 * pkin(4) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t41, t14, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t41, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t30, 0, -t27 * t21, -t30 * t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27, t30, 0, -t27 * pkin(8), -t30 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t9;
