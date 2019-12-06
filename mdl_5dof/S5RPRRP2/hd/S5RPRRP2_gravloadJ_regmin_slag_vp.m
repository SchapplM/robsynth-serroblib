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
% taug_reg [5x16]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:02
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
t11 = qJ(1) + pkin(8);
t12 = -qJ(5) - pkin(7);
t10 = qJ(3) + t11;
t7 = sin(t10);
t8 = cos(t10);
t15 = cos(qJ(4));
t9 = t15 * pkin(4) + pkin(3);
t20 = t7 * t12 - t8 * t9;
t4 = g(2) * t8 + g(3) * t7;
t3 = g(2) * t7 - g(3) * t8;
t19 = -t8 * t12 - t7 * t9;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t18 = g(2) * t16 + g(3) * t14;
t13 = sin(qJ(4));
t17 = -g(1) * t15 - t3 * t13;
t2 = t4 * t15;
t1 = t4 * t13;
t5 = [0, t18, -g(2) * t14 + g(3) * t16, t18 * pkin(1), 0, t4, -t3, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * (-pkin(2) * cos(t11) - t16 * pkin(1) + t20) - g(3) * (-pkin(2) * sin(t11) - t14 * pkin(1) + t19); 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1); 0, 0, 0, 0, 0, t4, -t3, 0, 0, 0, 0, 0, t2, -t1, t3, -g(2) * t20 - g(3) * t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, g(1) * t13 - t3 * t15, 0, t17 * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t4;];
taug_reg = t5;
