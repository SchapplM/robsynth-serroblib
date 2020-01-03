% Calculate minimal parameter regressor of joint inertia matrix for
% S5RRPRR3
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
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRPRR3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t23 = sin(pkin(9));
t37 = pkin(2) * t23;
t25 = sin(qJ(5));
t36 = pkin(4) * t25;
t20 = cos(qJ(2)) * pkin(1);
t19 = t20 + pkin(2);
t24 = cos(pkin(9));
t34 = sin(qJ(2)) * pkin(1);
t13 = t23 * t19 + t24 * t34;
t26 = sin(qJ(4));
t29 = cos(qJ(4));
t11 = t24 * t19 - t23 * t34;
t8 = pkin(3) + t11;
t6 = t29 * t8;
t4 = -t26 * t13 + t6;
t2 = -pkin(4) - t4;
t28 = cos(qJ(5));
t35 = t2 * t28;
t18 = t24 * pkin(2) + pkin(3);
t16 = t29 * t18;
t12 = -t26 * t37 + t16;
t9 = -pkin(4) - t12;
t33 = t9 * t28;
t32 = -t13 - t37;
t5 = -t29 * t13 - t26 * t8;
t14 = -t26 * t18 - t29 * t37;
t22 = t25 ^ 2;
t21 = pkin(4) * t28;
t17 = 0.2e1 * t25 * t28;
t10 = pkin(8) - t14;
t7 = t9 * t25;
t3 = pkin(8) - t5;
t1 = t2 * t25;
t15 = [1, 0, 0, 1, 0.2e1 * t20, -0.2e1 * t34, t11 ^ 2 + t13 ^ 2, 1, 0.2e1 * t4, 0.2e1 * t5, t22, t17, 0, 0, 0, -0.2e1 * t35, 0.2e1 * t1; 0, 0, 0, 1, t20, -t34, (t11 * t24 + t13 * t23) * pkin(2), 1, t32 * t26 + t16 + t6, t32 * t29 + (-t18 - t8) * t26, t22, t17, 0, 0, 0, (-t2 - t9) * t28, t7 + t1; 0, 0, 0, 1, 0, 0, (t23 ^ 2 + t24 ^ 2) * pkin(2) ^ 2, 1, 0.2e1 * t12, 0.2e1 * t14, t22, t17, 0, 0, 0, -0.2e1 * t33, 0.2e1 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, t4, t5, t22, t17, 0, 0, 0, t21 - t35, t1 - t36; 0, 0, 0, 0, 0, 0, 0, 1, t12, t14, t22, t17, 0, 0, 0, t21 - t33, t7 - t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t22, t17, 0, 0, 0, 0.2e1 * t21, -0.2e1 * t36; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * t3, -t28 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * t10, -t28 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t28, -t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t28, 0, -t25 * pkin(8), -t28 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg = t15;
