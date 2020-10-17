% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRPP4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:50
% EndTime: 2019-12-31 18:14:50
% DurationCPUTime: 0.16s
% Computational Cost: add. (106->43), mult. (157->45), div. (0->0), fcn. (143->6), ass. (0->23)
t21 = sin(qJ(1));
t23 = cos(qJ(1));
t34 = -g(1) * t21 + g(2) * t23;
t20 = sin(qJ(3));
t31 = t20 * pkin(3);
t30 = t23 * pkin(1) + t21 * qJ(2);
t15 = t23 * qJ(2);
t29 = -t21 * pkin(1) + t15;
t6 = g(1) * t23 + g(2) * t21;
t18 = qJ(3) + pkin(7);
t12 = sin(t18);
t13 = cos(t18);
t28 = t12 * pkin(4) - t13 * qJ(5);
t19 = -qJ(4) - pkin(6);
t27 = t21 * t19 + t23 * t31 + t29;
t26 = -t23 * t19 + t21 * t31 + t30;
t22 = cos(qJ(3));
t24 = g(3) * t20 + t22 * t34;
t4 = t6 * t13;
t3 = t6 * t12;
t2 = -g(3) * t12 - t34 * t13;
t1 = g(3) * t13 - t12 * t34;
t5 = [0, 0, 0, 0, 0, 0, -t34, t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, -t6, -g(1) * t29 - g(2) * t30, 0, 0, 0, 0, 0, 0, -t6 * t20, -t6 * t22, -t34, -g(1) * (t15 + (-pkin(1) - pkin(6)) * t21) - g(2) * (t23 * pkin(6) + t30), 0, 0, 0, 0, 0, 0, -t3, -t4, -t34, -g(1) * t27 - g(2) * t26, 0, 0, 0, 0, 0, 0, -t3, -t34, t4, -g(1) * (t28 * t23 + t27) - g(2) * (t28 * t21 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t24, g(3) * t22 - t20 * t34, 0, 0, 0, 0, 0, 0, 0, 0, -t2, t1, 0, t24 * pkin(3), 0, 0, 0, 0, 0, 0, -t2, 0, -t1, -g(3) * (-t28 - t31) + t34 * (pkin(3) * t22 + pkin(4) * t13 + qJ(5) * t12); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2;];
taug_reg = t5;
