% Calculate minimal parameter regressor of gravitation load for
% S5PRRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
% 
% Output:
% taug_reg [5x18]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5PRRRP4_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP4_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP4_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP4_gravloadJ_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
t24 = sin(qJ(4));
t26 = cos(qJ(4));
t46 = pkin(4) * t26 + qJ(5) * t24 + pkin(3);
t21 = qJ(2) + qJ(3);
t19 = sin(t21);
t22 = sin(pkin(8));
t23 = cos(pkin(8));
t32 = g(1) * t23 + g(2) * t22;
t30 = t32 * t19;
t25 = sin(qJ(2));
t45 = pkin(2) * t25;
t20 = cos(t21);
t43 = pkin(7) * t20;
t40 = g(3) * t19;
t39 = g(3) * t24;
t38 = t22 * t24;
t37 = t22 * t26;
t36 = t23 * t24;
t35 = t23 * t26;
t33 = t19 * pkin(7) + t46 * t20;
t6 = t20 * t38 + t35;
t8 = t20 * t36 - t37;
t1 = g(1) * t8 + g(2) * t6 + t19 * t39;
t7 = t20 * t37 - t36;
t9 = t20 * t35 + t38;
t29 = g(1) * t9 + g(2) * t7 + t26 * t40;
t4 = -g(3) * t20 + t30;
t28 = t46 * t30;
t27 = cos(qJ(2));
t12 = t23 * t43;
t10 = t22 * t43;
t5 = t32 * t20 + t40;
t3 = t4 * t26;
t2 = -t20 * t39 + t24 * t30;
t11 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, -g(3) * t27 + t32 * t25, g(3) * t25 + t32 * t27, 0, t4, t5, 0, 0, 0, 0, 0, t3, -t2, t3, -t5, t2, -g(1) * (-t23 * t45 + t12) - g(2) * (-t22 * t45 + t10) - g(3) * (t27 * pkin(2) + t33) + t28; 0, 0, 0, 0, 0, t4, t5, 0, 0, 0, 0, 0, t3, -t2, t3, -t5, t2, -g(1) * t12 - g(2) * t10 - g(3) * t33 + t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t29, t1, 0, -t29, -g(1) * (-t8 * pkin(4) + t9 * qJ(5)) - g(2) * (-t6 * pkin(4) + t7 * qJ(5)) - (-pkin(4) * t24 + qJ(5) * t26) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t11;
