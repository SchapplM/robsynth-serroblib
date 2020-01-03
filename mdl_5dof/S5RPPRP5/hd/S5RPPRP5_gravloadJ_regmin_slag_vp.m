% Calculate minimal parameter regressor of gravitation load for
% S5RPPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% taug_reg [5x22]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RPPRP5_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_gravloadJ_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = cos(pkin(7));
t25 = sin(qJ(4));
t23 = sin(pkin(7));
t37 = cos(qJ(4));
t32 = t23 * t37;
t11 = -t24 * t25 + t32;
t26 = sin(qJ(1));
t38 = g(1) * t26;
t27 = cos(qJ(1));
t35 = t24 * t27;
t34 = t27 * pkin(1) + t26 * qJ(2);
t33 = qJ(3) * t23;
t31 = pkin(2) * t35 + t27 * t33 + t34;
t4 = t11 * t26;
t6 = t25 * t35 - t27 * t32;
t30 = -g(1) * t4 - g(2) * t6;
t13 = g(1) * t27 + g(2) * t26;
t12 = -g(2) * t27 + t38;
t29 = -pkin(2) * t24 - pkin(1) - t33;
t10 = t23 * t25 + t24 * t37;
t1 = g(1) * t6 - g(2) * t4 + g(3) * t10;
t5 = t10 * t26;
t7 = t10 * t27;
t28 = g(1) * t7 + g(2) * t5 + g(3) * t11;
t20 = t27 * qJ(2);
t9 = t12 * t24;
t8 = t12 * t23;
t3 = g(3) * t24 - t13 * t23;
t2 = g(1) * t5 - g(2) * t7;
t14 = [0, t12, t13, t9, -t8, -t13, -g(1) * (-t26 * pkin(1) + t20) - g(2) * t34, t9, -t13, t8, -g(1) * t20 - g(2) * t31 - t29 * t38, 0, 0, 0, 0, 0, t2, -t30, t2, t13, t30, -g(1) * (-t5 * pkin(4) - t27 * pkin(6) + t4 * qJ(5) + t20) - g(2) * (pkin(3) * t35 + t7 * pkin(4) + t6 * qJ(5) + t31) + (-g(1) * (-pkin(3) * t24 + t29) + g(2) * pkin(6)) * t26; 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t28, t1, 0, -t28, -g(1) * (-t6 * pkin(4) + t7 * qJ(5)) - g(2) * (t4 * pkin(4) + t5 * qJ(5)) - g(3) * (-t10 * pkin(4) + t11 * qJ(5)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
