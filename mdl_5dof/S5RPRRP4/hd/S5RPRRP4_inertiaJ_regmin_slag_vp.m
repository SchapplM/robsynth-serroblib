% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 12:57
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP4_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 12:56:35
% EndTime: 2021-01-15 12:56:37
% DurationCPUTime: 0.32s
% Computational Cost: add. (360->67), mult. (761->120), div. (0->0), fcn. (806->6), ass. (0->45)
t30 = sin(qJ(4));
t32 = cos(qJ(4));
t28 = sin(pkin(8));
t33 = cos(qJ(3));
t42 = t33 * t28;
t31 = sin(qJ(3));
t44 = t31 * t28;
t12 = -t30 * t44 + t32 * t42;
t49 = -0.2e1 * t12;
t29 = cos(pkin(8));
t48 = 0.2e1 * t29;
t47 = t29 * pkin(4);
t46 = t30 * pkin(3);
t25 = t32 * pkin(3);
t16 = -t30 * t31 + t32 * t33;
t45 = t16 * t29;
t43 = t31 * t29;
t41 = t33 * t29;
t15 = pkin(3) * t44 + t28 * qJ(2);
t26 = t28 ^ 2;
t27 = t29 ^ 2;
t40 = t26 + t27;
t39 = qJ(2) * t31;
t38 = t26 * qJ(2);
t37 = qJ(2) * t41;
t18 = -t29 * pkin(2) - t28 * pkin(6) - pkin(1);
t14 = t33 * t18;
t6 = -pkin(7) * t42 + t14 + (-pkin(3) - t39) * t29;
t8 = t37 + (-pkin(7) * t28 + t18) * t31;
t3 = -t30 * t8 + t32 * t6;
t4 = t30 * t6 + t32 * t8;
t17 = t30 * t33 + t32 * t31;
t36 = -t12 * qJ(5) + t3;
t11 = t17 * t28;
t2 = -t11 * qJ(5) + t4;
t35 = 0.2e1 * pkin(4);
t23 = -0.2e1 * t46;
t22 = t25 + pkin(4);
t20 = t29 * t46;
t13 = t17 * t29;
t10 = t31 * t18 + t37;
t9 = -t29 * t39 + t14;
t7 = t11 * pkin(4) + t15;
t1 = t36 - t47;
t5 = [1, 0, 0, pkin(1) * t48, 0.2e1 * t40 * qJ(2), t40 * qJ(2) ^ 2 + pkin(1) ^ 2, t33 ^ 2 * t26, -0.2e1 * t33 * t26 * t31, -0.2e1 * t28 * t41, 0.2e1 * t28 * t43, t27, -0.2e1 * t9 * t29 + 0.2e1 * t31 * t38, 0.2e1 * t10 * t29 + 0.2e1 * t33 * t38, t12 ^ 2, t11 * t49, t29 * t49, t11 * t48, t27, 0.2e1 * t15 * t11 - 0.2e1 * t3 * t29, 0.2e1 * t15 * t12 + 0.2e1 * t4 * t29, -0.2e1 * t1 * t29 + 0.2e1 * t7 * t11, 0.2e1 * t7 * t12 + 0.2e1 * t2 * t29, -0.2e1 * t1 * t12 - 0.2e1 * t2 * t11, t1 ^ 2 + t2 ^ 2 + t7 ^ 2; 0, 0, 0, -t29, 0, -pkin(1), 0, 0, 0, 0, 0, -t41, t43, 0, 0, 0, 0, 0, -t45, t13, -t45, t13, -t17 * t11 - t16 * t12, t1 * t16 + t2 * t17; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16 ^ 2 + t17 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, t42, -t44, -t29, t9, -t10, 0, 0, t12, -t11, -t29, -t29 * t25 + t3, t20 - t4, (-pkin(4) - t22) * t29 + t36, -t2 + t20, -t11 * t46 - t22 * t12, t1 * t22 + t2 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t31, 0, 0, 0, 0, 0, t16, -t17, t16, -t17, 0, t16 * t22 + t17 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t25, t23, 0.2e1 * t22, t23, 0, t30 ^ 2 * pkin(3) ^ 2 + t22 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t11, -t29, t3, -t4, t36 - 0.2e1 * t47, -t2, -t12 * pkin(4), t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t17, t16, -t17, 0, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t25, -t46, t35 + t25, -t46, 0, t22 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t35, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, t12, 0, t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
