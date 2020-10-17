% Calculate minimal parameter regressor of joint inertia matrix for
% S4RRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% MM_reg [((4+1)*4/2)x21]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S4RRPR9_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR9_inertiaJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR9_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:10:04
% EndTime: 2019-12-31 17:10:05
% DurationCPUTime: 0.20s
% Computational Cost: add. (151->49), mult. (360->107), div. (0->0), fcn. (378->6), ass. (0->38)
t26 = sin(pkin(7));
t27 = cos(pkin(7));
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t13 = t28 * t26 - t30 * t27;
t29 = sin(qJ(2));
t10 = t13 * t29;
t43 = 0.2e1 * t10;
t21 = -t27 * pkin(3) - pkin(2);
t42 = 0.2e1 * t21;
t31 = cos(qJ(2));
t41 = 0.2e1 * t31;
t40 = pkin(5) * t26;
t22 = t29 * pkin(5);
t39 = t31 * pkin(5);
t38 = t26 * t29;
t37 = t27 * t29;
t36 = pkin(6) + qJ(3);
t16 = -t31 * pkin(2) - t29 * qJ(3) - pkin(1);
t8 = t26 * t16 + t27 * t39;
t35 = t26 ^ 2 + t27 ^ 2;
t12 = t27 * t16;
t7 = -t26 * t39 + t12;
t34 = -t7 * t26 + t8 * t27;
t33 = -pkin(2) * t29 + qJ(3) * t31;
t14 = t30 * t26 + t28 * t27;
t25 = t29 ^ 2;
t18 = t36 * t27;
t17 = t36 * t26;
t15 = pkin(3) * t38 + t22;
t9 = t14 * t29;
t6 = -t28 * t17 + t30 * t18;
t5 = -t30 * t17 - t28 * t18;
t4 = -pkin(6) * t38 + t8;
t3 = -pkin(6) * t37 + t12 + (-pkin(3) - t40) * t31;
t2 = t28 * t3 + t30 * t4;
t1 = -t28 * t4 + t30 * t3;
t11 = [1, 0, 0, t25, t29 * t41, 0, 0, 0, pkin(1) * t41, -0.2e1 * pkin(1) * t29, 0.2e1 * t25 * t40 - 0.2e1 * t7 * t31, 0.2e1 * t25 * pkin(5) * t27 + 0.2e1 * t8 * t31, 0.2e1 * (-t26 * t8 - t27 * t7) * t29, t25 * pkin(5) ^ 2 + t7 ^ 2 + t8 ^ 2, t10 ^ 2, t9 * t43, t31 * t43, t9 * t41, t31 ^ 2, -0.2e1 * t1 * t31 + 0.2e1 * t15 * t9, -0.2e1 * t15 * t10 + 0.2e1 * t2 * t31; 0, 0, 0, 0, 0, t29, t31, 0, -t22, -t39, -pkin(5) * t37 + t33 * t26, pkin(5) * t38 + t27 * t33, t34, -pkin(2) * t22 + qJ(3) * t34, -t10 * t14, t10 * t13 - t14 * t9, -t14 * t31, t13 * t31, 0, t15 * t13 + t21 * t9 - t5 * t31, -t21 * t10 + t15 * t14 + t6 * t31; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(2) * t27, -0.2e1 * pkin(2) * t26, 0.2e1 * t35 * qJ(3), qJ(3) ^ 2 * t35 + pkin(2) ^ 2, t14 ^ 2, -0.2e1 * t14 * t13, 0, 0, 0, t13 * t42, t14 * t42; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t37, 0, t22, 0, 0, 0, 0, 0, t9, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, t26, 0, -pkin(2), 0, 0, 0, 0, 0, t13, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, -t9, -t31, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, -t13, 0, t5, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t11;
