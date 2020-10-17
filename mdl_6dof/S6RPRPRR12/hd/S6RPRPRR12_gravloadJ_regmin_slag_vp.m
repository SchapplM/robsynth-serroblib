% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% taug_reg [6x31]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRR12_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR12_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:46:41
% EndTime: 2019-05-05 20:46:42
% DurationCPUTime: 0.21s
% Computational Cost: add. (116->48), mult. (214->70), div. (0->0), fcn. (227->8), ass. (0->35)
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t31 = t25 * pkin(3) - t28 * qJ(4);
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t45 = -g(1) * t26 + g(2) * t29;
t7 = g(3) * t28 - t25 * t45;
t42 = g(3) * t25;
t39 = t26 * t28;
t23 = qJ(5) + qJ(6);
t17 = sin(t23);
t38 = t29 * t17;
t18 = cos(t23);
t37 = t29 * t18;
t24 = sin(qJ(5));
t36 = t29 * t24;
t27 = cos(qJ(5));
t35 = t29 * t27;
t34 = t29 * pkin(1) + t26 * qJ(2);
t16 = g(1) * t29 + g(2) * t26;
t20 = t29 * qJ(2);
t14 = t16 * t28;
t13 = t16 * t25;
t12 = -t24 * t39 + t35;
t11 = -t27 * t39 - t36;
t10 = -t26 * t27 - t28 * t36;
t9 = t26 * t24 - t28 * t35;
t8 = -t28 * t45 - t42;
t6 = -t17 * t39 + t37;
t5 = -t18 * t39 - t38;
t4 = -t26 * t18 - t28 * t38;
t3 = t26 * t17 - t28 * t37;
t2 = g(1) * t6 - g(2) * t4 + t17 * t42;
t1 = -g(1) * t5 + g(2) * t3 - t18 * t42;
t15 = [0, -t45, t16, t45, -t16, -g(1) * (-t26 * pkin(1) + t20) - g(2) * t34, 0, 0, 0, 0, 0, -t13, -t14, -t45, t13, t14, -g(1) * (t31 * t29 + t20) - g(2) * (t29 * pkin(7) + t34) + (-g(1) * (-pkin(1) - pkin(7)) - g(2) * t31) * t26, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t12, -g(1) * t9 - g(2) * t11, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t45, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8, t7, 0, t8, -t7, g(3) * t31 + t45 * (pkin(3) * t28 + qJ(4) * t25) 0, 0, 0, 0, 0, -t7 * t24, -t7 * t27, 0, 0, 0, 0, 0, -t7 * t17, -t7 * t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t11 + g(2) * t9 - t27 * t42, g(1) * t12 - g(2) * t10 + t24 * t42, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t15;
