% Calculate minimal parameter regressor of gravitation load for
% S5RRPRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 21:03
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP10_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 21:02:45
% EndTime: 2021-01-15 21:02:46
% DurationCPUTime: 0.15s
% Computational Cost: add. (107->44), mult. (234->57), div. (0->0), fcn. (237->6), ass. (0->35)
t28 = sin(qJ(2));
t30 = cos(qJ(4));
t32 = cos(qJ(1));
t38 = t32 * t30;
t27 = sin(qJ(4));
t29 = sin(qJ(1));
t41 = t29 * t27;
t11 = t28 * t38 - t41;
t39 = t32 * t27;
t40 = t29 * t30;
t13 = t28 * t40 + t39;
t31 = cos(qJ(2));
t42 = g(3) * t31;
t1 = -g(1) * t11 - g(2) * t13 + t30 * t42;
t19 = g(1) * t32 + g(2) * t29;
t8 = g(3) * t28 + t19 * t31;
t37 = t31 * pkin(2) + t28 * qJ(3);
t22 = t27 * pkin(4) + qJ(3);
t26 = pkin(2) + pkin(7) + qJ(5);
t35 = t22 * t28 + t26 * t31;
t34 = g(1) * t29 - g(2) * t32;
t21 = t30 * pkin(4) + pkin(3) + pkin(6);
t17 = pkin(1) + t37;
t16 = t34 * t31;
t15 = t34 * t28;
t14 = -t28 * t41 + t38;
t12 = t28 * t39 + t40;
t9 = pkin(1) + t35;
t7 = t19 * t28 - t42;
t6 = t8 * t30;
t5 = t8 * t27;
t4 = -g(1) * t14 - g(2) * t12;
t3 = g(1) * t13 - g(2) * t11;
t2 = g(1) * t12 - g(2) * t14 - t27 * t42;
t10 = [0, t34, t19, 0, 0, 0, 0, 0, t16, -t15, -t19, -t16, t15, -g(1) * (t32 * pkin(6) - t17 * t29) - g(2) * (t29 * pkin(6) + t17 * t32), 0, 0, 0, 0, 0, t4, t3, t4, t3, t16, -g(1) * (t21 * t32 - t9 * t29) - g(2) * (t21 * t29 + t9 * t32); 0, 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, -t7, -t8, -g(3) * t37 - t19 * (-t28 * pkin(2) + t31 * qJ(3)), 0, 0, 0, 0, 0, -t5, -t6, -t5, -t6, t7, -g(3) * t35 - t19 * (t22 * t31 - t26 * t28); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t10;
