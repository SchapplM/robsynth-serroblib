% Calculate minimal parameter regressor of gravitation load for
% S5RRRRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
% 
% Output:
% taug_reg [5x23]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRR5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRR5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:22
% EndTime: 2020-01-03 12:13:23
% DurationCPUTime: 0.09s
% Computational Cost: add. (164->22), mult. (102->24), div. (0->0), fcn. (102->10), ass. (0->24)
t19 = qJ(1) + qJ(2);
t17 = qJ(3) + t19;
t11 = sin(t17);
t12 = cos(t17);
t24 = g(2) * t12 + g(3) * t11;
t7 = g(2) * t11 - g(3) * t12;
t23 = cos(qJ(1));
t22 = cos(qJ(4));
t21 = sin(qJ(1));
t20 = sin(qJ(4));
t18 = qJ(4) + qJ(5);
t16 = cos(t19);
t15 = cos(t18);
t14 = sin(t19);
t13 = sin(t18);
t10 = -g(2) * t16 - g(3) * t14;
t9 = g(2) * t14 - g(3) * t16;
t6 = t24 * t22;
t5 = t24 * t20;
t4 = t24 * t15;
t3 = t24 * t13;
t2 = -g(1) * t15 + t7 * t13;
t1 = g(1) * t13 + t7 * t15;
t8 = [0, -g(2) * t23 - g(3) * t21, g(2) * t21 - g(3) * t23, 0, t10, t9, 0, -t24, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, t10, t9, 0, -t24, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, -t24, t7, 0, 0, 0, 0, 0, -t6, t5, 0, 0, 0, 0, 0, -t4, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t22 + t7 * t20, g(1) * t20 + t7 * t22, 0, 0, 0, 0, 0, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1;];
taug_reg = t8;
