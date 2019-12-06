% Calculate minimal parameter regressor of gravitation load for
% S5RRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:31
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRPRR3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR3_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = qJ(1) + qJ(2);
t9 = pkin(9) + qJ(4) + t12;
t7 = sin(t9);
t8 = cos(t9);
t4 = g(2) * t8 + g(3) * t7;
t3 = -g(2) * t7 + g(3) * t8;
t10 = sin(t12);
t11 = cos(t12);
t6 = g(2) * t11 + g(3) * t10;
t16 = cos(qJ(1));
t15 = cos(qJ(5));
t14 = sin(qJ(1));
t13 = sin(qJ(5));
t5 = -g(2) * t10 + g(3) * t11;
t2 = t4 * t15;
t1 = t4 * t13;
t17 = [0, g(2) * t16 + g(3) * t14, -g(2) * t14 + g(3) * t16, 0, t6, t5, -g(2) * (-t16 * pkin(1) - pkin(2) * t11) - g(3) * (-t14 * pkin(1) - pkin(2) * t10), 0, t4, t3, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, t6, t5, t6 * pkin(2), 0, t4, t3, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t4, t3, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + t3 * t13, g(1) * t13 + t3 * t15;];
taug_reg = t17;
