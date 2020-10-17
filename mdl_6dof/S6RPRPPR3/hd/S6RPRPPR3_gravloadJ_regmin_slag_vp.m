% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:42:45
% EndTime: 2019-05-05 16:42:46
% DurationCPUTime: 0.20s
% Computational Cost: add. (168->53), mult. (199->70), div. (0->0), fcn. (195->8), ass. (0->37)
t25 = sin(qJ(3));
t23 = qJ(1) + pkin(9);
t17 = sin(t23);
t18 = cos(t23);
t9 = g(1) * t18 + g(2) * t17;
t48 = t9 * t25;
t19 = t25 * qJ(4);
t28 = cos(qJ(3));
t36 = t28 * pkin(3) + t19;
t2 = g(3) * t25 + t9 * t28;
t46 = pkin(3) * t25;
t45 = g(1) * t17;
t41 = g(3) * t28;
t40 = t28 * pkin(4);
t39 = t18 * t28;
t24 = sin(qJ(6));
t38 = t24 * t25;
t27 = cos(qJ(6));
t37 = t25 * t27;
t35 = qJ(4) * t28;
t26 = sin(qJ(1));
t34 = -t26 * pkin(1) + t18 * pkin(7);
t29 = cos(qJ(1));
t33 = t29 * pkin(1) + pkin(3) * t39 + t17 * pkin(7) + (pkin(2) + t19) * t18;
t32 = -g(2) * t18 + t45;
t31 = g(1) * t26 - g(2) * t29;
t30 = -pkin(2) - t36;
t12 = t18 * t35;
t10 = t17 * t35;
t8 = t32 * t28;
t7 = t32 * t25;
t6 = -t17 * t24 + t18 * t37;
t5 = -t17 * t27 - t18 * t38;
t4 = -t17 * t37 - t18 * t24;
t3 = t17 * t38 - t18 * t27;
t1 = -t41 + t48;
t11 = [0, t31, g(1) * t29 + g(2) * t26, t31 * pkin(1), 0, 0, 0, 0, 0, t8, -t7, t8, -t9, t7, -g(1) * t34 - g(2) * t33 - t30 * t45, t7, -t8, t9, -g(1) * (-t18 * qJ(5) + t34) - g(2) * (pkin(4) * t39 + t33) + (-g(1) * (t30 - t40) + g(2) * qJ(5)) * t17, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t6, -g(1) * t3 - g(2) * t5; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t18 * t46 + t12) - g(2) * (-t17 * t46 + t10) - g(3) * t36, -t2, -t1, 0, -g(1) * t12 - g(2) * t10 - g(3) * (t36 + t40) + (pkin(3) + pkin(4)) * t48, 0, 0, 0, 0, 0, -t2 * t27, t2 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t5 + g(2) * t3 - t24 * t41, g(1) * t6 - g(2) * t4 - t27 * t41;];
taug_reg  = t11;
