% Calculate inertial parameters regressor of gravitation load for
% S5PRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:38
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRPPR5_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR5_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRPPR5_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR5_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t23 = sin(qJ(2));
t20 = sin(pkin(7));
t21 = cos(pkin(7));
t30 = g(1) * t21 + g(2) * t20;
t43 = t30 * t23;
t25 = cos(qJ(2));
t34 = sin(pkin(8));
t35 = cos(pkin(8));
t7 = t23 * t34 + t25 * t35;
t27 = -t23 * t35 + t25 * t34;
t42 = g(3) * t27;
t40 = pkin(2) * t23;
t37 = t25 * pkin(2) + t23 * qJ(3);
t36 = qJ(3) * t25;
t33 = t25 * pkin(3) + t37;
t3 = t7 * t20;
t5 = t7 * t21;
t29 = -g(1) * t5 - g(2) * t3 + t42;
t4 = t27 * t20;
t6 = t27 * t21;
t28 = g(1) * t6 + g(2) * t4 + g(3) * t7;
t26 = (pkin(2) + pkin(3)) * t43;
t24 = cos(qJ(5));
t22 = sin(qJ(5));
t15 = t21 * t36;
t14 = t20 * t36;
t9 = g(1) * t20 - g(2) * t21;
t2 = g(3) * t23 + t30 * t25;
t1 = -g(3) * t25 + t43;
t8 = [0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t2, -g(1) * (-t21 * t40 + t15) - g(2) * (-t20 * t40 + t14) - g(3) * t37, 0, 0, 0, 0, 0, 0, -t28, t29, 0, -g(1) * t15 - g(2) * t14 - g(3) * t33 + t26, 0, 0, 0, 0, 0, 0, -t28 * t24, t28 * t22, -t29, -g(1) * (t6 * pkin(4) - t5 * pkin(6) + t15) - g(2) * (t4 * pkin(4) - t3 * pkin(6) + t14) - g(3) * (t7 * pkin(4) + pkin(6) * t27 + t33) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (-t20 * t24 - t5 * t22) - g(2) * (t21 * t24 - t3 * t22) - t22 * t42, -g(1) * (t20 * t22 - t5 * t24) - g(2) * (-t21 * t22 - t3 * t24) - t24 * t42, 0, 0;];
taug_reg = t8;
