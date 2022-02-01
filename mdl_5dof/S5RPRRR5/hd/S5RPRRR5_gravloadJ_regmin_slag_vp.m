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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:49:02
% EndTime: 2022-01-20 09:49:03
% DurationCPUTime: 0.08s
% Computational Cost: add. (120->15), mult. (78->21), div. (0->0), fcn. (76->8), ass. (0->20)
t11 = qJ(1) + pkin(9) + qJ(3);
t10 = cos(t11);
t9 = sin(t11);
t7 = g(1) * t9 - g(2) * t10;
t8 = g(1) * t10 + g(2) * t9;
t16 = sin(qJ(1));
t18 = cos(qJ(1));
t19 = g(1) * t16 - g(2) * t18;
t17 = cos(qJ(4));
t15 = sin(qJ(4));
t14 = qJ(4) + qJ(5);
t13 = cos(t14);
t12 = sin(t14);
t6 = t7 * t17;
t5 = t7 * t15;
t4 = t7 * t13;
t3 = t7 * t12;
t2 = g(3) * t12 + t8 * t13;
t1 = -g(3) * t13 + t8 * t12;
t20 = [0, t19, g(1) * t18 + g(2) * t16, t19 * pkin(1), 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t7, t8, 0, 0, 0, 0, 0, t6, -t5, 0, 0, 0, 0, 0, t4, -t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t17 + t8 * t15, g(3) * t15 + t8 * t17, 0, 0, 0, 0, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg = t20;
