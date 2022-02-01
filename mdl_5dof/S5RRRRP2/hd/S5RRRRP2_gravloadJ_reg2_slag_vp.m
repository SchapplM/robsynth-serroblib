% Calculate inertial parameters regressor of gravitation load for
% S5RRRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-20 11:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP2_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP2_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 11:49:28
% EndTime: 2022-01-20 11:49:28
% DurationCPUTime: 0.19s
% Computational Cost: add. (249->48), mult. (209->55), div. (0->0), fcn. (186->8), ass. (0->35)
t32 = -pkin(8) - pkin(7);
t29 = sin(qJ(1));
t42 = t29 * pkin(1);
t27 = qJ(1) + qJ(2);
t20 = sin(t27);
t22 = cos(t27);
t41 = t22 * pkin(2) + t20 * pkin(7);
t26 = qJ(3) + qJ(4);
t21 = cos(t26);
t30 = cos(qJ(3));
t23 = t30 * pkin(3);
t40 = pkin(4) * t21 + t23;
t39 = -t20 * pkin(2) + t22 * pkin(7);
t10 = pkin(2) + t40;
t25 = qJ(5) - t32;
t38 = t22 * t10 + t20 * t25;
t37 = -t20 * t10 + t25 * t22;
t18 = t23 + pkin(2);
t36 = t22 * t18 - t20 * t32;
t9 = g(1) * t22 + g(2) * t20;
t8 = g(1) * t20 - g(2) * t22;
t31 = cos(qJ(1));
t35 = g(1) * t29 - g(2) * t31;
t34 = -t20 * t18 - t22 * t32;
t19 = sin(t26);
t1 = -g(3) * t21 + t9 * t19;
t28 = sin(qJ(3));
t33 = -g(3) * t30 + t9 * t28;
t24 = t31 * pkin(1);
t6 = t8 * t30;
t5 = t8 * t28;
t4 = t8 * t21;
t3 = t8 * t19;
t2 = g(3) * t19 + t9 * t21;
t7 = [0, 0, 0, 0, 0, 0, t35, g(1) * t31 + g(2) * t29, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, t35 * pkin(1), 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * (t39 - t42) - g(2) * (t24 + t41), 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * (t34 - t42) - g(2) * (t24 + t36), 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * (t37 - t42) - g(2) * (t24 + t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t8, t9, 0, 0, 0, 0, 0, 0, 0, 0, t6, -t5, -t9, -g(1) * t39 - g(2) * t41, 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * t34 - g(2) * t36, 0, 0, 0, 0, 0, 0, t4, -t3, -t9, -g(1) * t37 - g(2) * t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, g(3) * t28 + t9 * t30, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t33 * pkin(3), 0, 0, 0, 0, 0, 0, t1, t2, 0, -g(3) * t40 - t9 * (-t28 * pkin(3) - pkin(4) * t19); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t8;];
taug_reg = t7;
