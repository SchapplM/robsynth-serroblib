% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x24]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 19:26
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP12_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 19:25:53
% EndTime: 2021-01-15 19:25:55
% DurationCPUTime: 0.32s
% Computational Cost: add. (198->65), mult. (369->119), div. (0->0), fcn. (352->4), ass. (0->45)
t44 = 2 * pkin(4);
t20 = cos(qJ(4));
t43 = 0.2e1 * t20;
t42 = 2 * qJ(2);
t18 = sin(qJ(4));
t41 = t18 * pkin(4);
t19 = sin(qJ(3));
t40 = t18 * t19;
t39 = t18 * t20;
t21 = cos(qJ(3));
t38 = t18 * t21;
t22 = -pkin(1) - pkin(6);
t37 = t18 * t22;
t36 = t19 * t22;
t35 = t20 * t19;
t12 = t20 * t21;
t34 = t20 * t22;
t33 = t21 * t19;
t32 = t21 * t22;
t31 = -qJ(5) - pkin(7);
t14 = t18 ^ 2;
t16 = t20 ^ 2;
t30 = t14 + t16;
t15 = t19 ^ 2;
t17 = t21 ^ 2;
t29 = -t15 - t17;
t28 = qJ(5) * t21;
t27 = -0.2e1 * t33;
t26 = t19 * t34;
t9 = t19 * pkin(3) - t21 * pkin(7) + qJ(2);
t5 = t20 * t9;
t25 = -t20 * t28 + t5;
t24 = -pkin(3) * t21 - pkin(7) * t19;
t10 = t31 * t18;
t11 = t31 * t20;
t23 = -t10 * t18 - t11 * t20;
t13 = -t20 * pkin(4) - pkin(3);
t8 = t29 * t20;
t7 = t29 * t18;
t6 = (-t22 + t41) * t21;
t4 = t18 * t9 + t26;
t3 = -t18 * t36 + t5;
t2 = t26 + (t9 - t28) * t18;
t1 = (pkin(4) - t37) * t19 + t25;
t45 = [1, 0, 0, -2 * pkin(1), t42, pkin(1) ^ 2 + qJ(2) ^ 2, t17, t27, 0, 0, 0, t19 * t42, t21 * t42, t16 * t17, -0.2e1 * t17 * t39, t33 * t43, t18 * t27, t15, -0.2e1 * t17 * t37 + 0.2e1 * t3 * t19, -0.2e1 * t17 * t34 - 0.2e1 * t4 * t19, 0.2e1 * t1 * t19 + 0.2e1 * t6 * t38, 0.2e1 * t6 * t12 - 0.2e1 * t2 * t19, 0.2e1 * (-t1 * t20 - t18 * t2) * t21, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, t7, t8, 0, -t6 * t21 + (-t1 * t18 + t2 * t20) * t19; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30 * t15 + t17; 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t32, -t36, t18 * t12, (-t14 + t16) * t21, t40, t35, 0, t24 * t18 + t20 * t32, -t18 * t32 + t24 * t20, t10 * t19 + t13 * t38 - t6 * t20, t11 * t19 + t13 * t12 + t6 * t18, (-t10 * t21 + t2) * t20 + (t11 * t21 - t1) * t18, t1 * t10 - t2 * t11 + t6 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, 0, 0, 0, 0, t12, -t38, t12, -t38, t30 * t19, -t13 * t21 + t23 * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t14, 0.2e1 * t39, 0, 0, 0, pkin(3) * t43, -0.2e1 * pkin(3) * t18, -0.2e1 * t13 * t20, 0.2e1 * t13 * t18, 0.2e1 * t23, t10 ^ 2 + t11 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, -t38, t19, t3, -t4, (t44 - t37) * t19 + t25, -t2, -pkin(4) * t12, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t40, -t35, -t40, -t35, 0, -pkin(4) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t18, t20, 0, -t18 * pkin(7), -t20 * pkin(7), t10, t11, -t41, t10 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t44, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t12, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t20, t18, 0, t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t45;
