% Calculate minimal parameter regressor of gravitation load for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPPR7_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPPR7_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 17:16:28
% EndTime: 2019-05-05 17:16:28
% DurationCPUTime: 0.19s
% Computational Cost: add. (117->42), mult. (179->61), div. (0->0), fcn. (175->8), ass. (0->31)
t21 = sin(qJ(1));
t24 = cos(qJ(1));
t41 = -g(1) * t21 + g(2) * t24;
t17 = qJ(3) + pkin(9);
t11 = sin(t17);
t38 = g(3) * t11;
t20 = sin(qJ(3));
t37 = t20 * pkin(3);
t19 = sin(qJ(6));
t36 = t21 * t19;
t22 = cos(qJ(6));
t35 = t21 * t22;
t34 = t24 * t19;
t33 = t24 * t22;
t32 = t24 * pkin(1) + t21 * qJ(2);
t31 = -t21 * pkin(1) + t24 * qJ(2);
t7 = g(1) * t24 + g(2) * t21;
t12 = cos(t17);
t30 = t11 * pkin(4) - t12 * qJ(5);
t18 = -qJ(4) - pkin(7);
t29 = t21 * t18 + t24 * t37 + t31;
t28 = -t24 * t18 + t21 * t37 + t32;
t26 = -g(3) * t12 + t11 * t41;
t23 = cos(qJ(3));
t25 = g(3) * t20 + t23 * t41;
t5 = -t12 * t36 + t33;
t4 = -t12 * t35 - t34;
t3 = -t12 * t34 - t35;
t2 = -t12 * t33 + t36;
t1 = -t12 * t41 - t38;
t6 = [0, -t41, t7, t41, -t7, -g(1) * t31 - g(2) * t32, 0, 0, 0, 0, 0, -t7 * t20, -t7 * t23, -t41, -g(1) * t29 - g(2) * t28, -t41, t7 * t11, t7 * t12, -g(1) * (t30 * t24 + t29) - g(2) * (t30 * t21 + t28) 0, 0, 0, 0, 0, -g(1) * t3 - g(2) * t5, -g(1) * t2 - g(2) * t4; 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, g(3) * t23 - t20 * t41, 0, t25 * pkin(3), 0, t1, t26, -g(3) * (-t30 - t37) + t41 * (pkin(3) * t23 + pkin(4) * t12 + qJ(5) * t11) 0, 0, 0, 0, 0, t26 * t19, t26 * t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t7, 0, 0, 0, -t7, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t4 + g(2) * t2 - t22 * t38, g(1) * t5 - g(2) * t3 + t19 * t38;];
taug_reg  = t6;
