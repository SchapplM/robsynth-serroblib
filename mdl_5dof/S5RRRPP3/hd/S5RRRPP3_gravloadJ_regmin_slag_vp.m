% Calculate minimal parameter regressor of gravitation load for
% S5RRRPP3
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
% taug_reg [5x21]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPP3_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP3_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP3_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP3_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t26 = sin(qJ(3));
t21 = t26 * qJ(4);
t47 = pkin(2) + t21;
t25 = qJ(1) + qJ(2);
t19 = sin(t25);
t20 = cos(t25);
t6 = g(1) * t20 + g(2) * t19;
t46 = t6 * t26;
t45 = g(1) * t19;
t28 = cos(qJ(3));
t22 = t28 * pkin(3);
t42 = t19 * t26;
t41 = t20 * t26;
t40 = t20 * t28;
t39 = pkin(3) + qJ(5);
t38 = t22 + t21;
t37 = qJ(4) * t28;
t36 = t28 * qJ(5);
t16 = t20 * pkin(7);
t27 = sin(qJ(1));
t35 = -t27 * pkin(1) + t16;
t33 = pkin(3) * t40 + t19 * pkin(7) + t47 * t20;
t32 = t19 * pkin(4) + t20 * t36 + t33;
t31 = (-t47 - t22) * t45;
t30 = (-t39 * t28 - t47) * t45;
t29 = cos(qJ(1));
t23 = t29 * pkin(1);
t17 = t20 * pkin(4);
t10 = t20 * t37;
t7 = t19 * t37;
t5 = -g(2) * t20 + t45;
t4 = -g(2) * t40 + t28 * t45;
t3 = g(1) * t42 - g(2) * t41;
t2 = g(3) * t26 + t6 * t28;
t1 = -g(3) * t28 + t46;
t8 = [0, g(1) * t27 - g(2) * t29, g(1) * t29 + g(2) * t27, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, -t6, -t4, t3, -g(1) * t35 - g(2) * (t23 + t33) - t31, -t6, t3, t4, -g(1) * (t17 + t35) - g(2) * (t23 + t32) - t30; 0, 0, 0, 0, t5, t6, 0, 0, 0, 0, 0, t4, -t3, -t6, -t4, t3, -g(1) * t16 - g(2) * t33 - t31, -t6, t3, t4, -g(1) * (t16 + t17) - g(2) * t32 - t30; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t2, 0, -t1, -t2, -g(1) * (-pkin(3) * t41 + t10) - g(2) * (-pkin(3) * t42 + t7) - g(3) * t38, 0, -t2, t1, -g(1) * t10 - g(2) * t7 - g(3) * (t36 + t38) + t39 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1, 0, 0, 0, -t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2;];
taug_reg = t8;
