% Calculate minimal parameter regressor of gravitation load for
% S5RPRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
% 
% Output:
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:24
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:24:01
% EndTime: 2021-01-15 11:24:02
% DurationCPUTime: 0.13s
% Computational Cost: add. (99->41), mult. (145->49), div. (0->0), fcn. (135->8), ass. (0->26)
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t32 = -g(1) * t21 + g(2) * t23;
t17 = qJ(3) + pkin(7);
t12 = sin(t17);
t29 = g(3) * t12;
t20 = sin(qJ(3));
t28 = t20 * pkin(3);
t16 = pkin(1) + pkin(6) + qJ(4);
t27 = t16 * t21;
t9 = g(1) * t23 + g(2) * t21;
t22 = cos(qJ(3));
t18 = sin(pkin(7));
t19 = cos(pkin(7));
t6 = pkin(4) * t19 + qJ(5) * t18 + pkin(3);
t7 = -t18 * pkin(4) + qJ(5) * t19;
t25 = t6 * t20 - t7 * t22 + qJ(2);
t24 = g(3) * t20 + t22 * t32;
t13 = cos(t17);
t11 = qJ(2) + t28;
t10 = t16 * t23;
t4 = t9 * t13;
t3 = t9 * t12;
t2 = t13 * t32 + t29;
t1 = g(3) * t13 - t12 * t32;
t5 = [0, -t32, t9, t32, -t9, -g(1) * (-t21 * pkin(1) + t23 * qJ(2)) - g(2) * (t23 * pkin(1) + t21 * qJ(2)), 0, 0, 0, 0, 0, -t9 * t20, -t9 * t22, -t3, -t4, -t32, -g(1) * (t11 * t23 - t27) - g(2) * (t11 * t21 + t10), -t3, -t32, t4, -g(1) * (t25 * t23 - t27) - g(2) * (t25 * t21 + t10); 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, g(3) * t22 - t20 * t32, t2, t1, 0, t24 * pkin(3), t2, 0, -t1, -g(3) * (-t12 * pkin(4) + t13 * qJ(5) - t28) + t32 * (t7 * t20 + t6 * t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t29 - t32 * (-t18 * t20 + t19 * t22);];
taug_reg = t5;
