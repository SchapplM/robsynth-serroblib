% Calculate minimal parameter regressor of gravitation load for
% S4RRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-01-31 13:16
% Revision: 9ef80adae39e3cd5824e7abdb6e4e1e7895c437e (2019-01-31)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4RRPR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR1_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t10 = qJ(1) + qJ(2);
t8 = sin(t10);
t9 = cos(t10);
t3 = g(1) * t8 - g(2) * t9;
t12 = cos(qJ(1));
t11 = sin(qJ(1));
t7 = pkin(7) + qJ(4) + t10;
t6 = cos(t7);
t5 = sin(t7);
t4 = g(1) * t9 + g(2) * t8;
t2 = g(1) * t6 + g(2) * t5;
t1 = g(1) * t5 - g(2) * t6;
t13 = [0, g(1) * t11 - g(2) * t12, g(1) * t12 + g(2) * t11, 0, t3, t4, -g(1) * (-pkin(1) * t11 - pkin(2) * t8) - g(2) * (pkin(1) * t12 + pkin(2) * t9) 0, t1, t2; 0, 0, 0, 0, t3, t4, t3 * pkin(2), 0, t1, t2; 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t13;
