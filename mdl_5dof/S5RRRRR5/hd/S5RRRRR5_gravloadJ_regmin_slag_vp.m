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
% Datum: 2022-01-20 12:02
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 12:02:06
% EndTime: 2022-01-20 12:02:07
% DurationCPUTime: 0.09s
% Computational Cost: add. (164->19), mult. (102->24), div. (0->0), fcn. (102->10), ass. (0->24)
t19 = qJ(1) + qJ(2);
t17 = qJ(3) + t19;
t11 = sin(t17);
t12 = cos(t17);
t8 = g(1) * t12 + g(2) * t11;
t7 = g(1) * t11 - g(2) * t12;
t23 = cos(qJ(1));
t22 = cos(qJ(4));
t21 = sin(qJ(1));
t20 = sin(qJ(4));
t18 = qJ(4) + qJ(5);
t16 = cos(t19);
t15 = cos(t18);
t14 = sin(t19);
t13 = sin(t18);
t10 = g(1) * t16 + g(2) * t14;
t9 = g(1) * t14 - g(2) * t16;
t6 = t7 * t22;
t5 = t7 * t20;
t4 = t7 * t15;
t3 = t7 * t13;
t2 = g(3) * t13 + t8 * t15;
t1 = -g(3) * t15 + t8 * t13;
t24 = [0, g(1) * t21 - g(2) * t23, g(1) * t23 + g(2) * t21, 0, t9, t10, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, t9, t10, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t22 + t8 * t20, g(3) * t20 + t8 * t22, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t24;
