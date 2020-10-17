% Calculate minimal parameter regressor of gravitation load for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP8_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:02:14
% EndTime: 2019-05-05 18:02:14
% DurationCPUTime: 0.22s
% Computational Cost: add. (183->59), mult. (280->76), div. (0->0), fcn. (291->8), ass. (0->38)
t25 = sin(qJ(1));
t28 = cos(qJ(1));
t51 = -g(1) * t25 + g(2) * t28;
t21 = qJ(3) + pkin(9);
t15 = sin(t21);
t16 = cos(t21);
t50 = -g(3) * t15 - t51 * t16;
t45 = g(3) * t16;
t44 = t16 * pkin(8);
t24 = sin(qJ(3));
t43 = t24 * pkin(3);
t23 = sin(qJ(5));
t42 = t25 * t23;
t26 = cos(qJ(5));
t41 = t25 * t26;
t40 = t28 * t23;
t39 = t28 * t26;
t38 = t28 * pkin(1) + t25 * qJ(2);
t37 = -t25 * pkin(1) + t28 * qJ(2);
t5 = t15 * t42 - t39;
t7 = t15 * t40 + t41;
t36 = g(1) * t7 + g(2) * t5;
t34 = t15 * pkin(4) - t44;
t10 = g(1) * t28 + g(2) * t25;
t22 = -qJ(4) - pkin(7);
t33 = t25 * t22 + t28 * t43 + t37;
t32 = -t28 * t22 + t25 * t43 + t38;
t31 = pkin(5) * t26 + qJ(6) * t23 + pkin(4);
t1 = g(1) * t5 - g(2) * t7 + t23 * t45;
t6 = t15 * t41 + t40;
t8 = t15 * t39 - t42;
t30 = g(1) * t6 - g(2) * t8 + t26 * t45;
t27 = cos(qJ(3));
t29 = g(3) * t24 + t27 * t51;
t4 = t50 * t26;
t3 = t50 * t23;
t2 = -g(1) * t8 - g(2) * t6;
t9 = [0, -t51, t10, t51, -t10, -g(1) * t37 - g(2) * t38, 0, 0, 0, 0, 0, -t10 * t24, -t10 * t27, -t51, -g(1) * t33 - g(2) * t32, 0, 0, 0, 0, 0, t2, t36, t2, t10 * t16, -t36, -g(1) * (t8 * pkin(5) + t7 * qJ(6) + t28 * t34 + t33) - g(2) * (t6 * pkin(5) + t5 * qJ(6) + t25 * t34 + t32); 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, t51, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, g(3) * t27 - t24 * t51, 0, t29 * pkin(3), 0, 0, 0, 0, 0, -t4, t3, -t4, t15 * t51 - t45, -t3, -g(3) * (-t15 * t31 - t43 + t44) + t51 * (pkin(3) * t27 + pkin(8) * t15 + t16 * t31); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t30, t1, 0, -t30, -g(1) * (-t5 * pkin(5) + t6 * qJ(6)) - g(2) * (t7 * pkin(5) - t8 * qJ(6)) - (-pkin(5) * t23 + qJ(6) * t26) * t45; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t9;
