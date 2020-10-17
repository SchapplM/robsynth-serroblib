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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:21
% DurationCPUTime: 0.26s
% Computational Cost: add. (148->52), mult. (277->99), div. (0->0), fcn. (260->4), ass. (0->41)
t18 = cos(qJ(4));
t40 = 0.2e1 * t18;
t39 = 2 * qJ(2);
t16 = sin(qJ(4));
t38 = t16 * pkin(4);
t17 = sin(qJ(3));
t37 = t16 * t17;
t36 = t16 * t18;
t19 = cos(qJ(3));
t35 = t16 * t19;
t20 = -pkin(1) - pkin(6);
t34 = t16 * t20;
t33 = t17 * t20;
t32 = t18 * t17;
t10 = t18 * t19;
t31 = t18 * t20;
t30 = t19 * t17;
t29 = t19 * t20;
t28 = -qJ(5) - pkin(7);
t12 = t16 ^ 2;
t14 = t18 ^ 2;
t27 = t12 + t14;
t13 = t17 ^ 2;
t15 = t19 ^ 2;
t26 = -t13 - t15;
t25 = qJ(5) * t19;
t24 = -0.2e1 * t30;
t23 = t17 * t31;
t22 = -pkin(3) * t19 - pkin(7) * t17;
t8 = t28 * t16;
t9 = t28 * t18;
t21 = -t8 * t16 - t9 * t18;
t11 = -t18 * pkin(4) - pkin(3);
t7 = t17 * pkin(3) - t19 * pkin(7) + qJ(2);
t6 = (-t20 + t38) * t19;
t5 = t18 * t7;
t4 = t16 * t7 + t23;
t3 = -t16 * t33 + t5;
t2 = t23 + (t7 - t25) * t16;
t1 = -t18 * t25 + t5 + (pkin(4) - t34) * t17;
t41 = [1, 0, 0, -2 * pkin(1), t39, pkin(1) ^ 2 + qJ(2) ^ 2, t15, t24, 0, 0, 0, t17 * t39, t19 * t39, t14 * t15, -0.2e1 * t15 * t36, t30 * t40, t16 * t24, t13, -0.2e1 * t15 * t34 + 0.2e1 * t3 * t17, -0.2e1 * t15 * t31 - 0.2e1 * t4 * t17, 0.2e1 * (-t1 * t18 - t16 * t2) * t19, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26 * t16, t26 * t18, 0, -t6 * t19 + (-t1 * t16 + t2 * t18) * t17; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t27 * t13 + t15; 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, t29, -t33, t16 * t10, (-t12 + t14) * t19, t37, t32, 0, t22 * t16 + t18 * t29, -t16 * t29 + t22 * t18, (-t19 * t8 + t2) * t18 + (t19 * t9 - t1) * t16, t1 * t8 + t6 * t11 - t2 * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -t17, 0, 0, 0, 0, 0, t10, -t35, t27 * t17, -t19 * t11 + t21 * t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t12, 0.2e1 * t36, 0, 0, 0, pkin(3) * t40, -0.2e1 * pkin(3) * t16, 0.2e1 * t21, t11 ^ 2 + t8 ^ 2 + t9 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t35, t17, t3, -t4, -pkin(4) * t10, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t37, -t32, 0, -pkin(4) * t37; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, t18, 0, -t16 * pkin(7), -t18 * pkin(7), -t38, t8 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t41;
