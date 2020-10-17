% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
% 
% Output:
% taug_reg [6x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:23:32
% EndTime: 2019-05-05 17:23:32
% DurationCPUTime: 0.18s
% Computational Cost: add. (88->60), mult. (211->73), div. (0->0), fcn. (207->6), ass. (0->35)
t25 = sin(qJ(3));
t28 = cos(qJ(3));
t40 = g(3) * t28;
t26 = sin(qJ(1));
t29 = cos(qJ(1));
t42 = g(2) * t29;
t9 = g(1) * t26 - t42;
t45 = t9 * t25 + t40;
t44 = -pkin(1) - pkin(7);
t43 = -pkin(3) - pkin(4);
t41 = g(3) * t25;
t39 = t25 * t26;
t38 = t25 * t29;
t37 = t26 * t28;
t24 = sin(qJ(6));
t36 = t29 * t24;
t27 = cos(qJ(6));
t35 = t29 * t27;
t32 = qJ(4) * t25;
t34 = pkin(3) * t37 + t26 * t32;
t33 = t29 * pkin(1) + t26 * qJ(2);
t18 = t28 * qJ(4);
t31 = pkin(3) * t39 + t29 * pkin(7) + t33;
t10 = g(1) * t29 + g(2) * t26;
t19 = t29 * qJ(2);
t30 = pkin(3) * t38 - t29 * t18 + t19;
t8 = t10 * t28;
t7 = t10 * t25;
t6 = t26 * t24 - t28 * t35;
t5 = t26 * t27 + t28 * t36;
t4 = t27 * t37 + t36;
t3 = t24 * t37 - t35;
t2 = g(1) * t37 - t28 * t42 - t41;
t1 = g(1) * t39 - g(2) * t38 + t40;
t11 = [0, t9, t10, -t9, -t10, -g(1) * (-t26 * pkin(1) + t19) - g(2) * t33, 0, 0, 0, 0, 0, -t7, -t8, -t7, t9, t8, -g(1) * (t44 * t26 + t30) - g(2) * (-t26 * t18 + t31) t8, t7, -t9, -g(1) * (pkin(4) * t38 + t30) - g(2) * (-t29 * qJ(5) + t31) + (-g(1) * (qJ(5) + t44) - g(2) * (t25 * pkin(4) - t18)) * t26, 0, 0, 0, 0, 0, -g(1) * t6 + g(2) * t4, -g(1) * t5 - g(2) * t3; 0, 0, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t9, 0, 0, 0, -t9, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, -t2, 0, -t1, -g(1) * t34 - g(3) * (-t25 * pkin(3) + t18) - (-pkin(3) * t28 - t32) * t42, -t1, t2, 0, -g(1) * (pkin(4) * t37 + t34) - g(3) * (t43 * t25 + t18) - (t43 * t28 - t32) * t42, 0, 0, 0, 0, 0, -t45 * t27, t45 * t24; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, 0, 0, t2, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t3 + g(2) * t5 + t24 * t41, -g(1) * t4 - g(2) * t6 + t27 * t41;];
taug_reg  = t11;
