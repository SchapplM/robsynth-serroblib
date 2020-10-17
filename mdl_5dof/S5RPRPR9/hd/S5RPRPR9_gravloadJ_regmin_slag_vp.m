% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:24:37
% EndTime: 2019-12-31 18:24:38
% DurationCPUTime: 0.12s
% Computational Cost: add. (104->35), mult. (132->49), div. (0->0), fcn. (133->8), ass. (0->27)
t12 = qJ(1) + pkin(8);
t10 = cos(t12);
t9 = sin(t12);
t25 = g(1) * t10 + g(2) * t9;
t14 = sin(qJ(3));
t17 = cos(qJ(3));
t2 = g(3) * t14 + t25 * t17;
t29 = g(3) * t17;
t13 = sin(qJ(5));
t28 = t13 * t14;
t16 = cos(qJ(5));
t27 = t14 * t16;
t26 = g(1) * t9 - g(2) * t10;
t15 = sin(qJ(1));
t18 = cos(qJ(1));
t24 = g(1) * t15 - g(2) * t18;
t23 = t17 * pkin(3) + t14 * qJ(4);
t21 = pkin(2) + t23;
t20 = t24 * pkin(1);
t8 = t26 * t17;
t7 = t26 * t14;
t6 = t10 * t16 - t9 * t28;
t5 = t10 * t13 + t9 * t27;
t4 = t10 * t28 + t9 * t16;
t3 = t10 * t27 - t9 * t13;
t1 = t25 * t14 - t29;
t11 = [0, t24, g(1) * t18 + g(2) * t15, t20, 0, 0, 0, 0, 0, t8, -t7, -t25, -t8, t7, t20 + (-g(2) * pkin(6) + g(1) * t21) * t9 + (-g(1) * pkin(6) - g(2) * t21) * t10, 0, 0, 0, 0, 0, -g(1) * t6 - g(2) * t4, g(1) * t5 - g(2) * t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(3) * t23 + t25 * (pkin(3) * t14 - qJ(4) * t17), 0, 0, 0, 0, 0, -t2 * t13, -t2 * t16; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5 + t16 * t29, g(1) * t4 - g(2) * t6 - t13 * t29;];
taug_reg = t11;
