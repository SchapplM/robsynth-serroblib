% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:00
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 22:59:35
% EndTime: 2021-01-15 22:59:36
% DurationCPUTime: 0.10s
% Computational Cost: add. (173->31), mult. (132->38), div. (0->0), fcn. (125->10), ass. (0->29)
t28 = cos(qJ(3));
t17 = t28 * pkin(3) + pkin(2);
t24 = qJ(1) + qJ(2);
t21 = sin(t24);
t22 = cos(t24);
t25 = -qJ(4) - pkin(7);
t32 = t21 * t17 + t22 * t25;
t23 = qJ(3) + pkin(9);
t31 = t22 * t17 - t21 * t25;
t10 = g(2) * t22 + g(3) * t21;
t9 = g(2) * t21 - g(3) * t22;
t26 = sin(qJ(3));
t30 = -g(1) * t28 + t9 * t26;
t29 = cos(qJ(1));
t27 = sin(qJ(1));
t20 = qJ(5) + t23;
t19 = cos(t23);
t18 = sin(t23);
t15 = cos(t20);
t14 = sin(t20);
t8 = t10 * t28;
t7 = t10 * t26;
t6 = t10 * t19;
t5 = t10 * t18;
t4 = t10 * t15;
t3 = t10 * t14;
t2 = -g(1) * t15 + t9 * t14;
t1 = g(1) * t14 + t9 * t15;
t11 = [0, -g(2) * t29 - g(3) * t27, g(2) * t27 - g(3) * t29, 0, -t10, t9, 0, 0, 0, 0, 0, -t8, t7, -t6, t5, -t9, -g(2) * (t29 * pkin(1) + t31) - g(3) * (t27 * pkin(1) + t32), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, -t8, t7, -t6, t5, -t9, -g(2) * t31 - g(3) * t32, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, g(1) * t26 + t9 * t28, -g(1) * t19 + t9 * t18, g(1) * t18 + t9 * t19, 0, t30 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t11;
