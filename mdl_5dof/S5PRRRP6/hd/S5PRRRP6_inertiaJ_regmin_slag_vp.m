% Calculate minimal parameter regressor of joint inertia matrix for
% S5PRRRP6
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
% MM_reg [((5+1)*5/2)x22]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5PRRRP6_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_inertiaJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_inertiaJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:52:17
% EndTime: 2019-12-05 16:52:19
% DurationCPUTime: 0.24s
% Computational Cost: add. (166->56), mult. (358->79), div. (0->0), fcn. (406->6), ass. (0->32)
t17 = sin(qJ(4));
t20 = cos(qJ(3));
t18 = sin(qJ(3));
t29 = cos(qJ(4));
t25 = t29 * t18;
t7 = t17 * t20 + t25;
t35 = -0.2e1 * t7;
t15 = -t20 * pkin(3) - pkin(2);
t34 = 0.2e1 * t15;
t33 = 0.2e1 * t20;
t32 = -pkin(7) - pkin(6);
t21 = cos(qJ(2));
t24 = t29 * t20;
t28 = t17 * t18;
t6 = -t24 + t28;
t31 = t21 * t6;
t30 = t21 * t7;
t19 = sin(qJ(2));
t27 = t18 * t19;
t26 = t29 * pkin(3);
t23 = 2 * pkin(4);
t22 = 2 * qJ(5);
t16 = t17 * pkin(3);
t13 = t26 + pkin(4);
t11 = t16 + qJ(5);
t9 = t32 * t20;
t5 = -t17 * t27 + t19 * t24;
t4 = t7 * t19;
t3 = t32 * t28 - t29 * t9;
t2 = -t17 * t9 - t32 * t25;
t1 = t6 * pkin(4) - t7 * qJ(5) + t15;
t8 = [1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21 ^ 2 + t4 ^ 2 + t5 ^ 2; 0, 0, t21, -t19, 0, 0, 0, 0, 0, t21 * t20, -t21 * t18, 0, 0, 0, 0, 0, -t31, -t30, -t31, t4 * t7 - t5 * t6, t30, -t21 * t1 + t4 * t2 + t5 * t3; 0, 1, 0, 0, t18 ^ 2, t18 * t33, 0, 0, 0, pkin(2) * t33, -0.2e1 * pkin(2) * t18, t7 ^ 2, t6 * t35, 0, 0, 0, t6 * t34, t7 * t34, 0.2e1 * t1 * t6, 0.2e1 * t2 * t7 - 0.2e1 * t3 * t6, t1 * t35, t1 ^ 2 + t2 ^ 2 + t3 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t27, -t20 * t19, 0, 0, 0, 0, 0, -t4, -t5, -t4, 0, t5, t5 * t11 - t4 * t13; 0, 0, 0, 0, 0, 0, t18, t20, 0, -t18 * pkin(6), -t20 * pkin(6), 0, 0, t7, -t6, 0, -t2, -t3, -t2, -t11 * t6 - t13 * t7, t3, t3 * t11 - t2 * t13; 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t26, -0.2e1 * t16, 0.2e1 * t13, 0, 0.2e1 * t11, t11 ^ 2 + t13 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4, -t5, -t4, 0, t5, -t4 * pkin(4) + t5 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, -t6, 0, -t2, -t3, -t2, -pkin(4) * t7 - t6 * qJ(5), t3, -t2 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t26, -t16, t23 + t26, 0, t22 + t16, t13 * pkin(4) + t11 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t23, 0, t22, pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t7, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t8;
