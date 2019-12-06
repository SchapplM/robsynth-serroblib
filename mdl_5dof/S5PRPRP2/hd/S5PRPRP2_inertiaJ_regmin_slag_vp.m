% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,theta1,theta3]';
% 
% Output:
% MM_reg [((5+1)*5/2)x17]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRPRP2_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRP2_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRP2_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t13 = sin(pkin(8));
t29 = -0.2e1 * t13;
t14 = cos(pkin(8));
t10 = t14 ^ 2;
t9 = t13 ^ 2;
t28 = t9 + t10;
t15 = sin(qJ(4));
t27 = t15 * t13;
t26 = t15 * t14;
t16 = cos(qJ(4));
t25 = t16 * t13;
t24 = t16 * t14;
t11 = t15 ^ 2;
t12 = t16 ^ 2;
t23 = t11 + t12;
t22 = qJ(3) * t15;
t21 = qJ(3) * t16;
t20 = qJ(5) * t13;
t19 = t14 * t21;
t7 = -t14 * pkin(3) - t13 * pkin(6) - pkin(2);
t5 = t16 * t7;
t1 = -t16 * t20 + t5 + (-pkin(4) - t22) * t14;
t2 = t19 + (t7 - t20) * t15;
t18 = t1 * t16 + t2 * t15;
t8 = t12 * t9;
t6 = (pkin(4) * t15 + qJ(3)) * t13;
t4 = t15 * t7 + t19;
t3 = -t14 * t22 + t5;
t17 = [1, 0, 0, 0, 0, 0, 0, t28, 0, 0, 0, 0, 0, 0, 0, 0, t11 * t9 + t10 + t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14 * t6 + (-t1 * t15 + t16 * t2) * t13; 0, 1, 0, 0, 0.2e1 * pkin(2) * t14, pkin(2) * t29, 0.2e1 * t28 * qJ(3), t28 * qJ(3) ^ 2 + pkin(2) ^ 2, t8, -0.2e1 * t16 * t9 * t15, t24 * t29, 0.2e1 * t13 * t26, t10, -0.2e1 * t3 * t14 + 0.2e1 * t9 * t22, 0.2e1 * t4 * t14 + 0.2e1 * t9 * t21, t18 * t29, t1 ^ 2 + t2 ^ 2 + t6 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, -t14, t13, 0, -pkin(2), 0, 0, 0, 0, 0, -t24, t26, -t23 * t13, t18; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t25, 0, -pkin(4) * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, -t27, -t14, t3, -t4, -pkin(4) * t25, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t16, -t15, 0, t16 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t17;
