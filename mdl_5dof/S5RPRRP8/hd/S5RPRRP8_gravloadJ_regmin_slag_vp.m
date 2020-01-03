% Calculate minimal parameter regressor of gravitation load for
% S5RPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% taug_reg [5x20]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:47
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPRRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP8_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP8_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t28 = sin(qJ(3));
t29 = sin(qJ(1));
t30 = cos(qJ(3));
t31 = cos(qJ(1));
t10 = t31 * t28 - t29 * t30;
t20 = sin(qJ(4));
t21 = cos(qJ(4));
t24 = -t21 * pkin(4) - t20 * qJ(5);
t22 = pkin(3) - t24;
t9 = -t29 * t28 - t31 * t30;
t40 = (g(2) * pkin(7) + g(1) * t22) * t10 - (-g(1) * pkin(7) + g(2) * t22) * t9;
t26 = g(1) * t10 - g(2) * t9;
t39 = t26 * t20;
t38 = t26 * t21;
t8 = g(1) * t9 + g(2) * t10;
t27 = t31 * pkin(1) + t29 * qJ(2);
t25 = -t29 * pkin(1) + t31 * qJ(2);
t12 = g(1) * t31 + g(2) * t29;
t11 = g(1) * t29 - g(2) * t31;
t2 = -g(3) * t20 - t8 * t21;
t1 = g(3) * t21 - t8 * t20;
t3 = [0, t11, t12, t11, -t12, -g(1) * t25 - g(2) * t27, 0, -t26, t8, 0, 0, 0, 0, 0, -t38, t39, -t38, -t8, -t39, -g(1) * (-t29 * pkin(2) + t25) - g(2) * (t31 * pkin(2) + t27) - t40; 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, t26, -t8, 0, 0, 0, 0, 0, t38, -t39, t38, t8, t39, t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(3) * t24 - t8 * (pkin(4) * t20 - qJ(5) * t21); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t3;
