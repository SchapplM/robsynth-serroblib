% Calculate minimal parameter regressor of gravitation load for
% S5RPRRR5
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:54:08
% EndTime: 2020-01-03 11:54:08
% DurationCPUTime: 0.08s
% Computational Cost: add. (120->17), mult. (78->21), div. (0->0), fcn. (76->8), ass. (0->20)
t11 = qJ(1) + pkin(9) + qJ(3);
t10 = cos(t11);
t9 = sin(t11);
t7 = g(2) * t9 - g(3) * t10;
t20 = g(2) * t10 + g(3) * t9;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t19 = -g(2) * t18 - g(3) * t16;
t17 = cos(qJ(4));
t15 = sin(qJ(4));
t14 = qJ(4) + qJ(5);
t13 = cos(t14);
t12 = sin(t14);
t6 = t20 * t17;
t5 = t20 * t15;
t4 = t20 * t13;
t3 = t20 * t12;
t2 = -g(1) * t13 + t7 * t12;
t1 = g(1) * t12 + t7 * t13;
t8 = [0, t19, g(2) * t16 - g(3) * t18, t19 * pkin(1), 0, -t20, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, -g(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, -t20, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t17 + t7 * t15, g(1) * t15 + t7 * t17, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t8;
