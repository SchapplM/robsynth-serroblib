% Calculate minimal parameter regressor of gravitation load for
% S5PRRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 15:53
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRPR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 15:52:41
% EndTime: 2021-01-15 15:52:42
% DurationCPUTime: 0.16s
% Computational Cost: add. (137->41), mult. (179->66), div. (0->0), fcn. (189->10), ass. (0->26)
t12 = sin(pkin(8));
t13 = cos(pkin(8));
t31 = -g(1) * t13 - g(2) * t12;
t16 = sin(qJ(2));
t18 = cos(qJ(2));
t3 = -g(3) * t18 - t31 * t16;
t28 = g(3) * t16;
t26 = t12 * t18;
t25 = t13 * t18;
t15 = sin(qJ(3));
t24 = t15 * t18;
t17 = cos(qJ(3));
t23 = t17 * t18;
t11 = qJ(3) + pkin(9);
t19 = -g(1) * (t12 * t17 - t13 * t24) - g(2) * (-t12 * t24 - t13 * t17) + t15 * t28;
t14 = qJ(4) + pkin(6);
t10 = qJ(5) + t11;
t9 = cos(t11);
t8 = sin(t11);
t7 = t17 * pkin(3) + pkin(2);
t6 = cos(t10);
t5 = sin(t10);
t4 = -t31 * t18 + t28;
t2 = -g(1) * (-t12 * t5 - t6 * t25) - g(2) * (t13 * t5 - t6 * t26) + t6 * t28;
t1 = -g(1) * (t12 * t6 - t5 * t25) - g(2) * (-t13 * t6 - t5 * t26) + t5 * t28;
t20 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, t3, t4, 0, 0, 0, 0, 0, t3 * t17, -t3 * t15, t3 * t9, -t3 * t8, -t4, -g(3) * (t16 * t14 + t18 * t7) + t31 * (t14 * t18 - t16 * t7), 0, 0, 0, 0, 0, t3 * t6, -t3 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, -g(1) * (-t12 * t15 - t13 * t23) - g(2) * (-t12 * t23 + t13 * t15) + t17 * t28, -g(1) * (t12 * t9 - t8 * t25) - g(2) * (-t13 * t9 - t8 * t26) + t8 * t28, -g(1) * (-t12 * t8 - t9 * t25) - g(2) * (t13 * t8 - t9 * t26) + t9 * t28, 0, t19 * pkin(3), 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t20;
