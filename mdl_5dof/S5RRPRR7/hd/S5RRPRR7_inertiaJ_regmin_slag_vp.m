% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% 
% Output:
% MM_reg [((5+1)*5/2)x23]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR7_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR7_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR7_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t38 = sin(qJ(2)) * pkin(1);
t17 = qJ(3) + t38;
t25 = sin(qJ(4));
t22 = t25 * pkin(4);
t13 = t17 + t22;
t42 = 0.2e1 * t13;
t41 = 0.2e1 * t17;
t19 = qJ(3) + t22;
t40 = 0.2e1 * t19;
t31 = 0.2e1 * qJ(3);
t24 = sin(qJ(5));
t39 = t24 * pkin(4);
t27 = cos(qJ(5));
t37 = t27 * pkin(4);
t28 = cos(qJ(4));
t36 = t28 * pkin(8);
t35 = cos(qJ(2)) * pkin(1);
t34 = t13 + t19;
t33 = qJ(3) + t17;
t20 = -pkin(2) - t35;
t32 = -0.2e1 * pkin(2);
t30 = -pkin(2) - pkin(7);
t23 = t28 ^ 2;
t21 = t28 * t30;
t16 = -0.2e1 * t28 * t25;
t15 = -pkin(7) + t20;
t14 = t28 * t15;
t12 = t21 - t36;
t11 = (-pkin(8) + t30) * t25;
t10 = -t24 * t25 + t27 * t28;
t9 = t24 * t28 + t27 * t25;
t8 = t10 ^ 2;
t7 = t14 - t36;
t6 = (-pkin(8) + t15) * t25;
t5 = -0.2e1 * t10 * t9;
t4 = -t27 * t11 - t24 * t12;
t3 = -t24 * t11 + t27 * t12;
t2 = -t24 * t7 - t27 * t6;
t1 = -t24 * t6 + t27 * t7;
t18 = [1, 0, 0, 1, 0.2e1 * t35, -0.2e1 * t38, 0.2e1 * t20, t41, t17 ^ 2 + t20 ^ 2, t23, t16, 0, 0, 0, t25 * t41, t28 * t41, t8, t5, 0, 0, 0, t9 * t42, t10 * t42; 0, 0, 0, 1, t35, -t38, t32 - t35, t31 + t38, -t20 * pkin(2) + t17 * qJ(3), t23, t16, 0, 0, 0, t33 * t25, t33 * t28, t8, t5, 0, 0, 0, t34 * t9, t34 * t10; 0, 0, 0, 1, 0, 0, t32, t31, pkin(2) ^ 2 + qJ(3) ^ 2, t23, t16, 0, 0, 0, t25 * t31, t28 * t31, t8, t5, 0, 0, 0, t9 * t40, t10 * t40; 0, 0, 0, 0, 0, 0, 1, 0, t20, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, -pkin(2), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, 0, t14, -t25 * t15, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, 0, t21, -t25 * t30, 0, 0, t10, -t9, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t37, -0.2e1 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, 0, t3, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t37, -t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t18;
