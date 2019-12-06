% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:36
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPPRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRP1_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t16 = sin(pkin(8));
t35 = -0.2e1 * t16;
t34 = 0.2e1 * t16;
t20 = sin(qJ(4));
t17 = sin(pkin(7));
t9 = t17 * pkin(1) + qJ(3);
t33 = t20 * t9;
t12 = t16 ^ 2;
t21 = cos(qJ(4));
t32 = t12 * t21;
t31 = t20 * t16;
t18 = cos(pkin(8));
t30 = t20 * t18;
t29 = t21 * t16;
t28 = t21 * t18;
t13 = t18 ^ 2;
t27 = t12 + t13;
t14 = t20 ^ 2;
t15 = t21 ^ 2;
t26 = t14 + t15;
t25 = qJ(5) * t16;
t24 = t9 * t28;
t19 = cos(pkin(7));
t11 = -t19 * pkin(1) - pkin(2);
t7 = -t18 * pkin(3) - t16 * pkin(6) + t11;
t5 = t21 * t7;
t1 = -t21 * t25 + t5 + (-pkin(4) - t33) * t18;
t2 = t24 + (t7 - t25) * t20;
t23 = t1 * t21 + t2 * t20;
t10 = t15 * t12;
t6 = (pkin(4) * t20 + t9) * t16;
t4 = t20 * t7 + t24;
t3 = -t9 * t30 + t5;
t8 = [1, 0, 0, (t17 ^ 2 + t19 ^ 2) * pkin(1) ^ 2, -0.2e1 * t11 * t18, t11 * t34, 0.2e1 * t27 * t9, t27 * t9 ^ 2 + t11 ^ 2, t10, -0.2e1 * t20 * t32, t28 * t35, t30 * t34, t13, 0.2e1 * t12 * t33 - 0.2e1 * t3 * t18, 0.2e1 * t4 * t18 + 0.2e1 * t9 * t32, t23 * t35, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6 * t18 + (-t1 * t20 + t2 * t21) * t16; 0, 0, 0, 1, 0, 0, 0, t27, 0, 0, 0, 0, 0, 0, 0, 0, t14 * t12 + t10 + t13; 0, 0, 0, 0, -t18, t16, 0, t11, 0, 0, 0, 0, 0, -t28, t30, -t26 * t16, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t31, -t18, t3, -t4, -pkin(4) * t29, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t29, 0, -pkin(4) * t31; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t20, 0, t21 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
