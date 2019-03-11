% Calculate inertial parameters regressor of gravitation load for
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
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRPRP8_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRP8_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_gravloadJ_reg2_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t34 = sin(qJ(1));
t37 = cos(qJ(1));
t63 = g(2) * t37;
t11 = g(1) * t34 - t63;
t30 = qJ(3) + pkin(9);
t25 = cos(t30);
t24 = sin(t30);
t62 = g(3) * t24;
t65 = t11 * t25 - t62;
t36 = cos(qJ(3));
t64 = pkin(3) * t36;
t61 = g(3) * t25;
t21 = t25 * pkin(8);
t33 = sin(qJ(3));
t60 = t33 * pkin(3);
t59 = t24 * t34;
t58 = t24 * t37;
t57 = t25 * t34;
t32 = sin(qJ(5));
t56 = t34 * t32;
t35 = cos(qJ(5));
t55 = t34 * t35;
t54 = t37 * t32;
t53 = t37 * t35;
t52 = pkin(1) * t37 + qJ(2) * t34;
t51 = pkin(4) * t57 + pkin(8) * t59 + t34 * t64;
t27 = t37 * qJ(2);
t50 = -pkin(1) * t34 + t27;
t49 = t21 - t60;
t7 = t24 * t56 - t53;
t9 = t24 * t54 + t55;
t48 = g(1) * t9 + g(2) * t7;
t47 = -pkin(8) * t24 - t64;
t12 = g(1) * t37 + g(2) * t34;
t46 = pkin(5) * t35 + qJ(6) * t32;
t31 = -qJ(4) - pkin(7);
t45 = t31 * t34 + t37 * t60 + t50;
t44 = -t31 * t37 + t34 * t60 + t52;
t43 = -pkin(4) - t46;
t1 = g(1) * t7 - g(2) * t9 + t32 * t61;
t10 = t24 * t53 - t56;
t8 = t24 * t55 + t54;
t42 = g(1) * t8 - g(2) * t10 + t35 * t61;
t40 = g(3) * t33 - t11 * t36;
t39 = pkin(4) * t58 - t21 * t37 + t45;
t38 = pkin(4) * t59 - pkin(8) * t57 + t44;
t6 = t12 * t25;
t5 = g(1) * t59 - g(2) * t58 + t61;
t4 = t65 * t35;
t3 = t65 * t32;
t2 = -g(1) * t10 - g(2) * t8;
t13 = [0, 0, 0, 0, 0, 0, t11, t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, -t12, -g(1) * t50 - g(2) * t52, 0, 0, 0, 0, 0, 0, -t12 * t33, -t12 * t36, t11, -g(1) * (t27 + (-pkin(1) - pkin(7)) * t34) - g(2) * (pkin(7) * t37 + t52) 0, 0, 0, 0, 0, 0, -t12 * t24, -t6, t11, -g(1) * t45 - g(2) * t44, 0, 0, 0, 0, 0, 0, t2, t48, t6, -g(1) * t39 - g(2) * t38, 0, 0, 0, 0, 0, 0, t2, t6, -t48, -g(1) * (pkin(5) * t10 + qJ(6) * t9 + t39) - g(2) * (pkin(5) * t8 + qJ(6) * t7 + t38); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t40, g(3) * t36 + t11 * t33, 0, 0, 0, 0, 0, 0, 0, 0, -t65, t5, 0, t40 * pkin(3), 0, 0, 0, 0, 0, 0, -t4, t3, -t5, -g(1) * t51 - g(3) * (-pkin(4) * t24 + t49) - (-pkin(4) * t25 + t47) * t63, 0, 0, 0, 0, 0, 0, -t4, -t5, -t3, -g(1) * (t46 * t57 + t51) - g(3) * t49 - t43 * t62 - (t25 * t43 + t47) * t63; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t12; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t42, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t42, -g(1) * (-pkin(5) * t7 + qJ(6) * t8) - g(2) * (pkin(5) * t9 - qJ(6) * t10) - (-pkin(5) * t32 + qJ(6) * t35) * t61; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t13;
