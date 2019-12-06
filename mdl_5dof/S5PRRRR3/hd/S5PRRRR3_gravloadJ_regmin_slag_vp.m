% Calculate minimal parameter regressor of gravitation load for
% S5PRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:06
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t15 = pkin(9) + qJ(2);
t14 = qJ(3) + t15;
t11 = qJ(4) + t14;
t7 = sin(t11);
t8 = cos(t11);
t4 = g(1) * t8 + g(2) * t7;
t3 = g(1) * t7 - g(2) * t8;
t17 = cos(qJ(5));
t16 = sin(qJ(5));
t13 = cos(t15);
t12 = sin(t15);
t10 = cos(t14);
t9 = sin(t14);
t6 = g(1) * t10 + g(2) * t9;
t5 = g(1) * t9 - g(2) * t10;
t2 = t3 * t17;
t1 = t3 * t16;
t18 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t12 - g(2) * t13, g(1) * t13 + g(2) * t12, 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t17 + t4 * t16, g(3) * t16 + t4 * t17;];
taug_reg = t18;
