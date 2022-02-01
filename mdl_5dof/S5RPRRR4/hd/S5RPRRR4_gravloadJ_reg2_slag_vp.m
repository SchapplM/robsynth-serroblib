% Calculate inertial parameters regressor of gravitation load for
% S5RPRRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR4_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:45
% DurationCPUTime: 0.14s
% Computational Cost: add. (238->38), mult. (116->44), div. (0->0), fcn. (102->10), ass. (0->29)
t21 = qJ(1) + pkin(9);
t19 = qJ(3) + t21;
t16 = qJ(4) + t19;
t11 = sin(t16);
t12 = cos(t16);
t33 = t12 * pkin(4) + t11 * pkin(8);
t14 = sin(t19);
t32 = pkin(3) * t14;
t18 = cos(t21);
t25 = cos(qJ(1));
t31 = t25 * pkin(1) + pkin(2) * t18;
t15 = cos(t19);
t10 = pkin(3) * t15;
t30 = t10 + t31;
t29 = -t11 * pkin(4) + t12 * pkin(8);
t17 = sin(t21);
t23 = sin(qJ(1));
t28 = -t23 * pkin(1) - pkin(2) * t17;
t4 = g(1) * t12 + g(2) * t11;
t3 = g(1) * t11 - g(2) * t12;
t5 = g(1) * t14 - g(2) * t15;
t27 = g(1) * t23 - g(2) * t25;
t26 = t28 - t32;
t24 = cos(qJ(5));
t22 = sin(qJ(5));
t6 = g(1) * t15 + g(2) * t14;
t2 = t3 * t24;
t1 = t3 * t22;
t7 = [0, 0, 0, 0, 0, 0, t27, g(1) * t25 + g(2) * t23, 0, 0, 0, 0, 0, 0, 0, 0, g(1) * t17 - g(2) * t18, g(1) * t18 + g(2) * t17, 0, t27 * pkin(1), 0, 0, 0, 0, 0, 0, t5, t6, 0, -g(1) * t28 - g(2) * t31, 0, 0, 0, 0, 0, 0, t3, t4, 0, -g(1) * t26 - g(2) * t30, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t26 + t29) - g(2) * (t30 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, t5 * pkin(3), 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * (t29 - t32) - g(2) * (t10 + t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, 0, 0, 0, t2, -t1, -t4, -g(1) * t29 - g(2) * t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t24 + t4 * t22, g(3) * t22 + t4 * t24, 0, 0;];
taug_reg = t7;
