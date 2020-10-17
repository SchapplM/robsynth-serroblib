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
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:09:29
% EndTime: 2020-01-03 12:09:29
% DurationCPUTime: 0.09s
% Computational Cost: add. (139->26), mult. (106->32), div. (0->0), fcn. (99->8), ass. (0->24)
t23 = cos(qJ(3));
t15 = t23 * pkin(3) + pkin(2);
t19 = qJ(1) + qJ(2);
t17 = sin(t19);
t18 = cos(t19);
t20 = -qJ(4) - pkin(7);
t27 = t17 * t15 + t18 * t20;
t26 = t18 * t15 - t17 * t20;
t8 = g(2) * t18 + g(3) * t17;
t7 = g(2) * t17 - g(3) * t18;
t21 = sin(qJ(3));
t25 = -g(1) * t23 + t7 * t21;
t24 = cos(qJ(1));
t22 = sin(qJ(1));
t16 = qJ(3) + pkin(9) + qJ(5);
t13 = cos(t16);
t12 = sin(t16);
t6 = t8 * t23;
t5 = t8 * t21;
t4 = t8 * t13;
t3 = t8 * t12;
t2 = -g(1) * t13 + t7 * t12;
t1 = g(1) * t12 + t7 * t13;
t9 = [0, -g(2) * t24 - g(3) * t22, g(2) * t22 - g(3) * t24, 0, -t8, t7, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * (t24 * pkin(1) + t26) - g(3) * (t22 * pkin(1) + t27), 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, -t8, t7, 0, 0, 0, 0, 0, -t6, t5, -t7, -g(2) * t26 - g(3) * t27, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(1) * t21 + t7 * t23, 0, t25 * pkin(3), 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t9;
