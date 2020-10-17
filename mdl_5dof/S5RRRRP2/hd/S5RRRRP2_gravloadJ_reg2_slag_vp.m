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
% Datum: 2020-01-03 12:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 12:11:46
% EndTime: 2020-01-03 12:11:46
% DurationCPUTime: 0.18s
% Computational Cost: add. (249->49), mult. (209->55), div. (0->0), fcn. (186->8), ass. (0->35)
t37 = -pkin(8) - pkin(7);
t31 = qJ(3) + qJ(4);
t25 = cos(t31);
t35 = cos(qJ(3));
t28 = t35 * pkin(3);
t43 = pkin(4) * t25 + t28;
t11 = pkin(2) + t43;
t32 = qJ(1) + qJ(2);
t24 = sin(t32);
t26 = cos(t32);
t30 = -qJ(5) + t37;
t46 = t24 * t11 + t26 * t30;
t22 = t28 + pkin(2);
t45 = t24 * t22 + t26 * t37;
t44 = t26 * pkin(2) + t24 * pkin(7);
t42 = t24 * pkin(2) - t26 * pkin(7);
t41 = t26 * t11 - t24 * t30;
t40 = t26 * t22 - t24 * t37;
t10 = g(2) * t26 + g(3) * t24;
t9 = g(2) * t24 - g(3) * t26;
t34 = sin(qJ(1));
t36 = cos(qJ(1));
t39 = -g(2) * t36 - g(3) * t34;
t23 = sin(t31);
t2 = -g(1) * t25 + t9 * t23;
t33 = sin(qJ(3));
t38 = -g(1) * t35 + t9 * t33;
t29 = t36 * pkin(1);
t27 = t34 * pkin(1);
t6 = t10 * t35;
t5 = t10 * t33;
t4 = t10 * t25;
t3 = t10 * t23;
t1 = g(1) * t23 + t9 * t25;
t7 = [0, 0, 0, 0, 0, 0, t39, g(2) * t34 - g(3) * t36, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, t39 * pkin(1), 0, 0, 0, 0, 0, 0, -t6, t5, -t9, -g(2) * (t29 + t44) - g(3) * (t27 + t42), 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * (t29 + t40) - g(3) * (t27 + t45), 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * (t29 + t41) - g(3) * (t27 + t46); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, t9, 0, 0, 0, 0, 0, 0, 0, 0, -t6, t5, -t9, -g(2) * t44 - g(3) * t42, 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * t40 - g(3) * t45, 0, 0, 0, 0, 0, 0, -t4, t3, -t9, -g(2) * t41 - g(3) * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, g(1) * t33 + t9 * t35, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t38 * pkin(3), 0, 0, 0, 0, 0, 0, t2, t1, 0, -g(1) * t43 - t9 * (-t33 * pkin(3) - pkin(4) * t23); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, 0, t2 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10;];
taug_reg = t7;
