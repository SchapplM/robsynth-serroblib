% Calculate minimal parameter regressor of joint inertia matrix for
% S5PPRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
% 
% Output:
% MM_reg [((5+1)*5/2)x16]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 15:11
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PPRRP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP3_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP3_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t11 = sin(qJ(4));
t30 = -0.2e1 * t11;
t13 = cos(qJ(4));
t29 = 0.2e1 * t13;
t6 = t11 ^ 2;
t28 = t13 ^ 2 + t6;
t27 = t11 * pkin(6);
t12 = sin(qJ(3));
t9 = sin(pkin(8));
t26 = t12 * t9;
t25 = t13 * pkin(6);
t14 = cos(qJ(3));
t24 = t14 * t9;
t23 = t11 * t12;
t22 = t13 * t12;
t21 = t14 * t11;
t4 = t14 * t13;
t20 = t9 * t23;
t19 = t9 * t22;
t18 = t28 * t12;
t10 = cos(pkin(8));
t1 = t10 * t13 + t9 * t21;
t2 = -t10 * t11 + t9 * t4;
t17 = t1 * t11 + t2 * t13;
t16 = -t11 * pkin(4) + t13 * qJ(5);
t7 = t12 ^ 2;
t5 = t9 ^ 2;
t3 = -t13 * pkin(4) - t11 * qJ(5) - pkin(3);
t8 = [1, t10 ^ 2 + t5, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1 ^ 2 + t2 ^ 2 + t7 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, (t17 - t24) * t12; 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14 ^ 2 + t28 * t7; 0, 0, 0, -t26, -t24, 0, 0, 0, 0, 0, -t19, t20, -t19, t17, -t20, t17 * pkin(6) + t3 * t26; 0, 0, 0, t14, -t12, 0, 0, 0, 0, 0, t4, -t21, t4, t18, t21, pkin(6) * t18 - t14 * t3; 0, 0, 1, 0, 0, t6, t11 * t29, 0, 0, 0, pkin(3) * t29, pkin(3) * t30, -0.2e1 * t3 * t13, 0.2e1 * t28 * pkin(6), t3 * t30, t28 * pkin(6) ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, -t2, -t1, 0, t2, -t1 * pkin(4) + t2 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23, -t22, -t23, 0, t22, t16 * t12; 0, 0, 0, 0, 0, 0, 0, t11, t13, 0, -t27, -t25, -t27, t16, t25, t16 * pkin(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t11, 0, t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
