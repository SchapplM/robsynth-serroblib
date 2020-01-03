% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% taug_reg [5x25]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t23 = cos(qJ(2));
t20 = qJ(2) + qJ(3);
t17 = sin(t20);
t18 = cos(t20);
t33 = t18 * pkin(3) + t17 * qJ(4);
t30 = t23 * pkin(2) + t33;
t41 = pkin(1) + t30;
t24 = cos(qJ(1));
t40 = g(2) * t24;
t22 = sin(qJ(1));
t35 = g(1) * t24;
t6 = g(2) * t22 + t35;
t39 = t6 * t17;
t21 = sin(qJ(2));
t37 = pkin(2) * t21;
t36 = pkin(3) * t17;
t13 = t18 * pkin(4);
t32 = qJ(4) * t18;
t25 = -pkin(7) - pkin(6);
t31 = qJ(5) + t25;
t29 = t41 * t40;
t28 = -t36 - t37;
t5 = g(1) * t22 - t40;
t26 = (pkin(3) + pkin(4)) * t39;
t9 = t24 * t32;
t7 = t22 * t32;
t4 = t5 * t18;
t3 = t5 * t17;
t2 = g(3) * t17 + t6 * t18;
t1 = -g(3) * t18 + t39;
t8 = [0, t5, t6, 0, 0, 0, 0, 0, t5 * t23, -t5 * t21, 0, 0, 0, 0, 0, t4, -t3, t4, -t6, t3, t25 * t35 - t29 + (g(1) * t41 + g(2) * t25) * t22, t4, t3, t6, -t29 + (g(1) * t31 - g(2) * t13) * t24 + (-g(1) * (-t41 - t13) + g(2) * t31) * t22; 0, 0, 0, 0, 0, 0, 0, 0, -g(3) * t23 + t6 * t21, g(3) * t21 + t6 * t23, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (t28 * t24 + t9) - g(2) * (t28 * t22 + t7) - g(3) * t30, t1, -t2, 0, -g(1) * (-t24 * t37 + t9) - g(2) * (-t22 * t37 + t7) - g(3) * (t13 + t30) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, t1, 0, -t2, -g(1) * (-t24 * t36 + t9) - g(2) * (-t22 * t36 + t7) - g(3) * t33, t1, -t2, 0, -g(1) * t9 - g(2) * t7 - g(3) * (t13 + t33) + t26; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t5;];
taug_reg = t8;
