% Calculate minimal parameter regressor of gravitation load for
% S5RPRPR15
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x24]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPR15_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:37:24
% EndTime: 2019-12-31 18:37:25
% DurationCPUTime: 0.18s
% Computational Cost: add. (86->41), mult. (163->63), div. (0->0), fcn. (168->8), ass. (0->31)
t19 = sin(qJ(3));
t21 = cos(qJ(3));
t24 = t19 * pkin(3) - t21 * qJ(4);
t20 = sin(qJ(1));
t22 = cos(qJ(1));
t41 = -g(1) * t20 + g(2) * t22;
t6 = -g(3) * t19 - t21 * t41;
t37 = g(3) * t21;
t16 = pkin(8) + qJ(5);
t10 = sin(t16);
t35 = t20 * t10;
t11 = cos(t16);
t34 = t20 * t11;
t17 = sin(pkin(8));
t33 = t20 * t17;
t18 = cos(pkin(8));
t32 = t20 * t18;
t31 = t22 * t10;
t30 = t22 * t11;
t29 = t22 * t17;
t28 = t22 * t18;
t27 = t22 * pkin(1) + t20 * qJ(2);
t9 = g(1) * t22 + g(2) * t20;
t13 = t22 * qJ(2);
t7 = t9 * t21;
t5 = -t19 * t41 + t37;
t4 = t19 * t30 - t35;
t3 = t19 * t31 + t34;
t2 = t19 * t34 + t31;
t1 = -t19 * t35 + t30;
t8 = [0, -t41, t9, t41, -t9, -g(1) * (-t20 * pkin(1) + t13) - g(2) * t27, 0, 0, 0, 0, 0, -t9 * t19, -t7, -g(1) * (t19 * t28 - t33) - g(2) * (t19 * t32 + t29), -g(1) * (-t19 * t29 - t32) - g(2) * (-t19 * t33 + t28), t7, -g(1) * (t24 * t22 + t13) - g(2) * (t22 * pkin(6) + t27) + (-g(1) * (-pkin(1) - pkin(6)) - g(2) * t24) * t20, 0, 0, 0, 0, 0, -g(1) * t4 - g(2) * t2, g(1) * t3 - g(2) * t1; 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t41, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t6 * t18, t6 * t17, -t5, g(3) * t24 + t41 * (pkin(3) * t21 + qJ(4) * t19), 0, 0, 0, 0, 0, -t6 * t11, t6 * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t6, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t3 + t10 * t37, g(1) * t2 - g(2) * t4 + t11 * t37;];
taug_reg = t8;
