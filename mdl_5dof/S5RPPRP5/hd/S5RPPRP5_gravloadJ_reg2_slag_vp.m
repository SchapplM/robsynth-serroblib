% Calculate inertial parameters regressor of gravitation load for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_gravloadJ_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:44
% DurationCPUTime: 0.17s
% Computational Cost: add. (108->48), mult. (265->64), div. (0->0), fcn. (287->6), ass. (0->33)
t25 = cos(pkin(7));
t26 = sin(qJ(4));
t24 = sin(pkin(7));
t41 = cos(qJ(4));
t36 = t24 * t41;
t11 = -t25 * t26 + t36;
t27 = sin(qJ(1));
t42 = g(1) * t27;
t28 = cos(qJ(1));
t39 = t25 * t28;
t38 = t28 * pkin(1) + t27 * qJ(2);
t37 = qJ(3) * t24;
t21 = t28 * qJ(2);
t35 = -t28 * pkin(6) + t21;
t34 = pkin(2) * t39 + t28 * t37 + t38;
t4 = t11 * t27;
t6 = t26 * t39 - t28 * t36;
t33 = -g(1) * t4 - g(2) * t6;
t32 = pkin(3) * t39 + t34;
t13 = g(1) * t28 + g(2) * t27;
t12 = -g(2) * t28 + t42;
t31 = -pkin(2) * t25 - pkin(1) - t37;
t10 = t24 * t26 + t25 * t41;
t1 = g(1) * t6 - g(2) * t4 + g(3) * t10;
t5 = t10 * t27;
t7 = t10 * t28;
t30 = g(1) * t7 + g(2) * t5 + g(3) * t11;
t29 = (-g(1) * (-pkin(3) * t25 + t31) + g(2) * pkin(6)) * t27;
t9 = t12 * t25;
t8 = t12 * t24;
t3 = g(3) * t25 - t13 * t24;
t2 = g(1) * t5 - g(2) * t7;
t14 = [0, 0, 0, 0, 0, 0, t12, t13, 0, 0, 0, 0, 0, 0, 0, 0, t9, -t8, -t13, -g(1) * (-t27 * pkin(1) + t21) - g(2) * t38, 0, 0, 0, 0, 0, 0, t9, -t13, t8, -g(1) * t21 - g(2) * t34 - t31 * t42, 0, 0, 0, 0, 0, 0, t2, -t33, t13, -g(1) * t35 - g(2) * t32 + t29, 0, 0, 0, 0, 0, 0, t2, t13, t33, -g(1) * (-t5 * pkin(4) + t4 * qJ(5) + t35) - g(2) * (t7 * pkin(4) + t6 * qJ(5) + t32) + t29; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t30, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t30, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t10 * pkin(4) + t11 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
