% Calculate minimal parameter regressor of gravitation load for
% S4PRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [4x10]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S4PRRR1_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRR1_gravloadJ_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4PRRR1_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4PRRR1_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t13 = pkin(7) + qJ(2);
t12 = qJ(3) + t13;
t11 = cos(t13);
t10 = sin(t13);
t9 = qJ(4) + t12;
t8 = cos(t12);
t7 = sin(t12);
t6 = cos(t9);
t5 = sin(t9);
t4 = g(1) * t8 + g(2) * t7;
t3 = g(1) * t7 - g(2) * t8;
t2 = g(1) * t6 + g(2) * t5;
t1 = g(1) * t5 - g(2) * t6;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, g(1) * t10 - g(2) * t11, g(1) * t11 + g(2) * t10, 0, t3, t4, 0, t1, t2; 0, 0, 0, 0, 0, t3, t4, 0, t1, t2; 0, 0, 0, 0, 0, 0, 0, 0, t1, t2;];
taug_reg  = t14;
