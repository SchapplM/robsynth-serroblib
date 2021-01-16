% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 16:46
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP7_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 16:45:23
% EndTime: 2021-01-15 16:45:25
% DurationCPUTime: 0.32s
% Computational Cost: add. (214->73), mult. (512->140), div. (0->0), fcn. (569->8), ass. (0->45)
t26 = sin(qJ(3));
t47 = -0.2e1 * t26;
t29 = cos(qJ(3));
t46 = 0.2e1 * t29;
t28 = cos(qJ(4));
t45 = pkin(3) * t28;
t25 = sin(qJ(4));
t44 = pkin(7) * t25;
t43 = t25 * pkin(4);
t24 = cos(pkin(5));
t23 = sin(pkin(5));
t41 = t23 * sin(qJ(2));
t10 = -t24 * t29 + t26 * t41;
t42 = t10 * t28;
t40 = t23 * cos(qJ(2));
t39 = t25 * t26;
t38 = t25 * t28;
t37 = t25 * t29;
t18 = t28 * t26;
t36 = t28 * t29;
t35 = -qJ(5) - pkin(8);
t34 = qJ(5) * t26;
t33 = t26 * t46;
t32 = pkin(7) * t36;
t14 = -t29 * pkin(3) - t26 * pkin(8) - pkin(2);
t12 = t28 * t14;
t31 = -t28 * t34 + t12;
t22 = t28 ^ 2;
t21 = t26 ^ 2;
t20 = t25 ^ 2;
t19 = -t28 * pkin(4) - pkin(3);
t16 = t35 * t28;
t15 = t35 * t25;
t13 = (pkin(7) + t43) * t26;
t11 = t24 * t26 + t29 * t41;
t9 = t10 * t25;
t8 = t25 * t14 + t32;
t7 = -pkin(7) * t37 + t12;
t6 = t32 + (t14 - t34) * t25;
t5 = t11 * t28 - t25 * t40;
t4 = -t11 * t25 - t28 * t40;
t3 = (-pkin(4) - t44) * t29 + t31;
t2 = t10 * t18 + t5 * t29;
t1 = t10 * t39 - t4 * t29;
t17 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, t40, -t41, 0, 0, 0, 0, 0, t29 * t40, -t26 * t40, 0, 0, 0, 0, 0, t1, t2, t1, t2, (-t25 * t5 - t28 * t4) * t26, t10 * t13 + t4 * t3 + t5 * t6; 0, 1, 0, 0, t21, t33, 0, 0, 0, pkin(2) * t46, pkin(2) * t47, t22 * t21, -0.2e1 * t21 * t38, t36 * t47, t25 * t33, t29 ^ 2, 0.2e1 * t21 * t44 - 0.2e1 * t7 * t29, 0.2e1 * t21 * pkin(7) * t28 + 0.2e1 * t8 * t29, 0.2e1 * t13 * t39 - 0.2e1 * t3 * t29, 0.2e1 * t13 * t18 + 0.2e1 * t6 * t29, 0.2e1 * (-t25 * t6 - t28 * t3) * t26, t13 ^ 2 + t3 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t11, 0, 0, 0, 0, 0, -t42, t9, -t42, t9, -t4 * t25 + t5 * t28, t10 * t19 + t4 * t15 - t5 * t16; 0, 0, 0, 0, 0, 0, t26, t29, 0, -t26 * pkin(7), -t29 * pkin(7), t25 * t18, (-t20 + t22) * t26, -t37, -t36, 0, -pkin(7) * t18 + (-pkin(3) * t26 + pkin(8) * t29) * t25, pkin(8) * t36 + (t44 - t45) * t26, -t13 * t28 - t15 * t29 + t19 * t39, t13 * t25 - t16 * t29 + t19 * t18, (-t15 * t26 + t6) * t28 + (t16 * t26 - t3) * t25, t13 * t19 + t3 * t15 - t6 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t20, 0.2e1 * t38, 0, 0, 0, 0.2e1 * t45, -0.2e1 * pkin(3) * t25, -0.2e1 * t19 * t28, 0.2e1 * t19 * t25, -0.2e1 * t15 * t25 - 0.2e1 * t16 * t28, t15 ^ 2 + t16 ^ 2 + t19 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, -t5, t4, -t5, 0, t4 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, -t39, -t29, t7, -t8, (-0.2e1 * pkin(4) - t44) * t29 + t31, -t6, -pkin(4) * t18, t3 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * pkin(8), -t28 * pkin(8), t15, t16, -t43, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39, t18, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t28, t25, 0, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t17;
