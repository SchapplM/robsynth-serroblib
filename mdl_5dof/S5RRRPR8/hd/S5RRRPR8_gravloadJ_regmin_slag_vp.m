% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:21
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR8_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:20:26
% EndTime: 2019-12-31 21:20:26
% DurationCPUTime: 0.18s
% Computational Cost: add. (150->46), mult. (200->60), div. (0->0), fcn. (200->8), ass. (0->36)
t26 = cos(qJ(2));
t21 = qJ(2) + qJ(3);
t18 = sin(t21);
t19 = cos(t21);
t34 = t19 * pkin(3) + t18 * qJ(4);
t42 = t26 * pkin(2) + t34;
t24 = sin(qJ(1));
t27 = cos(qJ(1));
t11 = g(1) * t27 + g(2) * t24;
t4 = g(3) * t18 + t11 * t19;
t41 = pkin(3) * t18;
t39 = g(3) * t19;
t22 = sin(qJ(5));
t38 = t24 * t22;
t25 = cos(qJ(5));
t37 = t24 * t25;
t36 = t27 * t22;
t35 = t27 * t25;
t33 = qJ(4) * t19;
t23 = sin(qJ(2));
t32 = -pkin(2) * t23 - t41;
t31 = g(1) * t24 - g(2) * t27;
t30 = pkin(1) + t42;
t28 = -pkin(7) - pkin(6);
t13 = t27 * t33;
t12 = t24 * t33;
t10 = -t18 * t38 + t35;
t9 = t18 * t37 + t36;
t8 = t18 * t36 + t37;
t7 = t18 * t35 - t38;
t6 = t31 * t19;
t5 = t31 * t18;
t3 = t11 * t18 - t39;
t2 = t4 * t25;
t1 = t4 * t22;
t14 = [0, t31, t11, 0, 0, 0, 0, 0, t31 * t26, -t31 * t23, 0, 0, 0, 0, 0, t6, -t5, -t11, -t6, t5, (g(1) * t28 - g(2) * t30) * t27 + (g(1) * t30 + g(2) * t28) * t24, 0, 0, 0, 0, 0, -g(1) * t10 - g(2) * t8, g(1) * t9 - g(2) * t7; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t26 + t11 * t23, g(3) * t23 + t11 * t26, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (t32 * t27 + t13) - g(2) * (t32 * t24 + t12) - g(3) * t42, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, -t3, -t4, -g(1) * (-t27 * t41 + t13) - g(2) * (-t24 * t41 + t12) - g(3) * t34, 0, 0, 0, 0, 0, -t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t3, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t7 - g(2) * t9 + t25 * t39, g(1) * t8 - g(2) * t10 - t22 * t39;];
taug_reg = t14;
