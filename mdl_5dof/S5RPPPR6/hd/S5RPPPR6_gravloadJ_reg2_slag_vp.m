% Calculate inertial parameters regressor of gravitation load for
% S5RPPPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta4]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPPR6_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR6_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPPR6_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR6_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t29 = sin(qJ(1));
t49 = g(1) * t29;
t27 = cos(pkin(7));
t28 = sin(qJ(5));
t48 = t27 * t28;
t30 = cos(qJ(5));
t47 = t27 * t30;
t31 = cos(qJ(1));
t46 = t27 * t31;
t24 = sin(pkin(8));
t45 = t29 * t24;
t26 = cos(pkin(8));
t44 = t29 * t26;
t43 = t31 * t24;
t42 = t31 * t26;
t19 = t31 * qJ(2);
t41 = t31 * pkin(3) + t19;
t40 = t31 * pkin(1) + t29 * qJ(2);
t25 = sin(pkin(7));
t39 = qJ(3) * t25;
t38 = -pkin(1) - t39;
t37 = pkin(2) * t46 + t31 * t39 + t40;
t5 = -t25 * t42 + t45;
t7 = t25 * t44 + t43;
t36 = g(1) * t7 + g(2) * t5;
t14 = g(1) * t31 + g(2) * t29;
t13 = -g(2) * t31 + t49;
t35 = t29 * pkin(3) + qJ(4) * t46 + t37;
t8 = -t25 * t45 + t42;
t34 = t8 * t28 + t29 * t47;
t33 = t29 * t48 - t8 * t30;
t32 = ((-pkin(2) - qJ(4)) * t27 + t38) * t49;
t10 = t13 * t27;
t9 = t13 * t25;
t6 = t25 * t43 + t44;
t4 = -g(3) * t25 - t14 * t27;
t3 = g(3) * t27 - t14 * t25;
t2 = t28 * t46 + t6 * t30;
t1 = -t6 * t28 + t30 * t46;
t11 = [0, 0, 0, 0, 0, 0, t13, t14, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t9, -t14, -g(1) * (-t29 * pkin(1) + t19) - g(2) * t40, 0, 0, 0, 0, 0, 0, -t14, -t10, t9, -g(1) * t19 - g(2) * t37 - (-pkin(2) * t27 + t38) * t49, 0, 0, 0, 0, 0, 0, -g(1) * t8 - g(2) * t6, t36, t10, -g(1) * t41 - g(2) * t35 - t32, 0, 0, 0, 0, 0, 0, g(1) * t33 - g(2) * t2, g(1) * t34 - g(2) * t1, -t36, -g(1) * (t8 * pkin(4) + t7 * pkin(6) + t41) - g(2) * (t6 * pkin(4) + t5 * pkin(6) + t35) - t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t13; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t1 - g(2) * t34 - g(3) * (t24 * t48 + t25 * t30), g(1) * t2 + g(2) * t33 - g(3) * (t24 * t47 - t25 * t28), 0, 0;];
taug_reg = t11;
