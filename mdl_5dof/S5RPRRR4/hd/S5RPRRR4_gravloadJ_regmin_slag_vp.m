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
% Datum: 2022-01-23 09:35
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% OptimizationMode: 2
% StartTime: 2022-01-23 09:34:45
% EndTime: 2022-01-23 09:34:46
% DurationCPUTime: 0.07s
% Computational Cost: add. (124->14), mult. (62->19), div. (0->0), fcn. (60->8), ass. (0->18)
t12 = qJ(1) + pkin(9) + qJ(3);
t11 = qJ(4) + t12;
t7 = sin(t11);
t8 = cos(t11);
t4 = g(1) * t8 + g(2) * t7;
t3 = g(1) * t7 - g(2) * t8;
t14 = sin(qJ(1));
t16 = cos(qJ(1));
t17 = g(1) * t14 - g(2) * t16;
t15 = cos(qJ(5));
t13 = sin(qJ(5));
t10 = cos(t12);
t9 = sin(t12);
t6 = g(1) * t10 + g(2) * t9;
t5 = g(1) * t9 - g(2) * t10;
t2 = t3 * t15;
t1 = t3 * t13;
t18 = [0, t17, g(1) * t16 + g(2) * t14, t17 * pkin(1), 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, t5, t6, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, t3, t4, 0, 0, 0, 0, 0, t2, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t15 + t4 * t13, g(3) * t13 + t4 * t15;];
taug_reg = t18;
