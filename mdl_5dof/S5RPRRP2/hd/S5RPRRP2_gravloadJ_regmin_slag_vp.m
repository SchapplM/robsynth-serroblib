% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2022-01-23 09:28
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP2_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:28:12
% EndTime: 2022-01-23 09:28:13
% DurationCPUTime: 0.09s
% Computational Cost: add. (146->27), mult. (102->29), div. (0->0), fcn. (93->8), ass. (0->20)
t12 = qJ(1) + pkin(8);
t16 = cos(qJ(4));
t10 = t16 * pkin(4) + pkin(3);
t13 = -qJ(5) - pkin(7);
t11 = qJ(3) + t12;
t8 = sin(t11);
t9 = cos(t11);
t20 = t9 * t10 - t8 * t13;
t6 = g(1) * t9 + g(2) * t8;
t5 = g(1) * t8 - g(2) * t9;
t15 = sin(qJ(1));
t17 = cos(qJ(1));
t19 = g(1) * t15 - g(2) * t17;
t18 = -t8 * t10 - t9 * t13;
t14 = sin(qJ(4));
t1 = -g(3) * t16 + t6 * t14;
t4 = t5 * t16;
t3 = t5 * t14;
t2 = g(3) * t14 + t6 * t16;
t7 = [0, t19, g(1) * t17 + g(2) * t15, t19 * pkin(1), 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * (-pkin(2) * sin(t12) - t15 * pkin(1) + t18) - g(2) * (pkin(2) * cos(t12) + t17 * pkin(1) + t20); 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, t4, -t3, -t6, -g(1) * t18 - g(2) * t20; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, t2, 0, t1 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t5;];
taug_reg = t7;
