% Calculate minimal parameter regressor of joint inertia matrix for
% S5RPRRP1
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
% Datum: 2019-12-05 18:00
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RPRRP1_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP1_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP1_inertiaJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
t18 = sin(qJ(4));
t19 = sin(qJ(3));
t20 = cos(qJ(4));
t21 = cos(qJ(3));
t9 = -t18 * t19 + t20 * t21;
t7 = t9 ^ 2;
t8 = t18 * t21 + t20 * t19;
t30 = t8 ^ 2 + t7;
t12 = t19 * pkin(3) + qJ(2);
t28 = 0.2e1 * t12;
t27 = 0.2e1 * qJ(2);
t26 = t9 * pkin(4);
t25 = t18 * pkin(3);
t22 = -pkin(1) - pkin(6);
t10 = (-pkin(7) + t22) * t19;
t15 = t21 * t22;
t11 = -t21 * pkin(7) + t15;
t3 = -t18 * t10 + t20 * t11;
t1 = -t9 * qJ(5) + t3;
t4 = -t20 * t10 - t18 * t11;
t2 = -t8 * qJ(5) - t4;
t24 = t1 * t9 + t2 * t8;
t17 = t20 * pkin(3);
t14 = t17 + pkin(4);
t23 = t9 * t14 + t8 * t25;
t5 = t8 * pkin(4) + t12;
t6 = [1, 0, 0, -2 * pkin(1), t27, (pkin(1) ^ 2) + qJ(2) ^ 2, t21 ^ 2, -0.2e1 * t21 * t19, 0, 0, 0, t19 * t27, t21 * t27, t7, -0.2e1 * t9 * t8, 0, 0, 0, t8 * t28, t9 * t28, -0.2e1 * t24, t1 ^ 2 + t2 ^ 2 + t5 ^ 2; 0, 0, 0, 1, 0, -pkin(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t24; 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30; 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, t15, -t19 * t22, 0, 0, t9, -t8, 0, t3, t4, -t23, t1 * t14 + t2 * t25; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, -t19, 0, 0, 0, 0, 0, t9, -t8, 0, t23; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t17, -0.2e1 * t25, 0, t18 ^ 2 * pkin(3) ^ 2 + t14 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t3, t4, -t26, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, 0, t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t17, -t25, 0, t14 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, pkin(4) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t6;
