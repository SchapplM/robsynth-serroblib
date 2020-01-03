% Calculate minimal parameter regressor of gravitation load for
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
% taug_reg [5x17]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR4_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t12 = qJ(1) + pkin(9) + qJ(3);
t11 = qJ(4) + t12;
t7 = sin(t11);
t8 = cos(t11);
t18 = g(2) * t8 + g(3) * t7;
t3 = g(2) * t7 - g(3) * t8;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t17 = -g(2) * t16 - g(3) * t14;
t15 = cos(qJ(5));
t13 = sin(qJ(5));
t10 = cos(t12);
t9 = sin(t12);
t6 = -g(2) * t10 - g(3) * t9;
t5 = g(2) * t9 - g(3) * t10;
t2 = t18 * t15;
t1 = t18 * t13;
t4 = [0, t17, g(2) * t14 - g(3) * t16, t17 * pkin(1), 0, t6, t5, 0, -t18, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t6, t5, 0, -t18, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, -t18, t3, 0, 0, 0, 0, 0, -t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t15 + t3 * t13, g(1) * t13 + t3 * t15;];
taug_reg = t4;
