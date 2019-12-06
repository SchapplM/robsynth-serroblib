% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% MM_reg [((5+1)*5/2)x20]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t22 = cos(qJ(3));
t16 = -t22 * pkin(3) - pkin(2);
t29 = 0.2e1 * t16;
t28 = 0.2e1 * t22;
t27 = pkin(6) + pkin(7);
t18 = sin(qJ(4));
t26 = t18 * pkin(3);
t19 = sin(qJ(3));
t20 = sin(qJ(2));
t25 = t19 * t20;
t24 = t22 * t20;
t11 = t27 * t19;
t12 = t27 * t22;
t21 = cos(qJ(4));
t3 = -t21 * t11 - t18 * t12;
t4 = t18 * t11 - t21 * t12;
t9 = t18 * t22 + t21 * t19;
t23 = cos(qJ(2));
t17 = t21 * pkin(3);
t15 = t17 + pkin(4);
t8 = t18 * t19 - t21 * t22;
t7 = -t18 * t25 + t21 * t24;
t6 = t9 * t20;
t5 = t8 * pkin(4) + t16;
t2 = -t8 * qJ(5) - t4;
t1 = -t9 * qJ(5) + t3;
t10 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t23 ^ 2 + t6 ^ 2 + t7 ^ 2; 0, 0, t23, -t20, 0, 0, 0, 0, 0, t23 * t22, -t23 * t19, 0, 0, 0, 0, 0, -t23 * t8, -t23 * t9, t6 * t9 - t7 * t8, -t6 * t1 + t7 * t2 - t23 * t5; 0, 1, 0, 0, t19 ^ 2, t19 * t28, 0, 0, 0, pkin(2) * t28, -0.2e1 * pkin(2) * t19, t9 ^ 2, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t29, t9 * t29, -0.2e1 * t1 * t9 - 0.2e1 * t2 * t8, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t25, -t24, 0, 0, 0, 0, 0, -t6, -t7, 0, -t6 * t15 + t7 * t26; 0, 0, 0, 0, 0, 0, t19, t22, 0, -t19 * pkin(6), -t22 * pkin(6), 0, 0, t9, -t8, 0, t3, t4, -t15 * t9 - t8 * t26, t1 * t15 + t2 * t26; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t26, 0, t18 ^ 2 * pkin(3) ^ 2 + t15 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, -t7, 0, -t6 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t3, t4, -pkin(4) * t9, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t26, 0, t15 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t10;
