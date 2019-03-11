% Calculate inertial parameters regressor of gravitation load for
% S6RPRRRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 06:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRRP3_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP3_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP3_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP3_gravloadJ_reg2_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
t43 = cos(qJ(4));
t32 = t43 * pkin(4) + pkin(3);
t41 = sin(qJ(3));
t44 = cos(qJ(3));
t46 = -pkin(9) - pkin(8);
t58 = t44 * t32 - t41 * t46;
t38 = qJ(1) + pkin(10);
t33 = sin(t38);
t34 = cos(t38);
t20 = g(1) * t34 + g(2) * t33;
t49 = -g(3) * t44 + t20 * t41;
t80 = g(1) * t33;
t77 = g(3) * t41;
t40 = sin(qJ(4));
t75 = t33 * t40;
t74 = t33 * t43;
t73 = t34 * t40;
t72 = t34 * t43;
t39 = qJ(4) + qJ(5);
t35 = sin(t39);
t71 = t35 * t41;
t70 = t35 * t44;
t36 = cos(t39);
t69 = t36 * t41;
t68 = t36 * t44;
t67 = t40 * t44;
t65 = t43 * t44;
t64 = t44 * t46;
t63 = t34 * t67;
t45 = cos(qJ(1));
t62 = t45 * pkin(1) + t34 * pkin(2) + t33 * pkin(7);
t42 = sin(qJ(1));
t61 = -t42 * pkin(1) + t34 * pkin(7);
t10 = t33 * t68 - t34 * t35;
t9 = t33 * t70 + t34 * t36;
t60 = -t9 * pkin(5) + t10 * qJ(6);
t11 = -t33 * t36 + t34 * t70;
t12 = t33 * t35 + t34 * t68;
t59 = -t11 * pkin(5) + t12 * qJ(6);
t57 = g(1) * t9 - g(2) * t11;
t56 = t44 * pkin(3) + t41 * pkin(8);
t54 = -g(2) * t34 + t80;
t53 = g(1) * t42 - g(2) * t45;
t52 = pkin(5) * t36 + qJ(6) * t35;
t14 = t33 * t67 + t72;
t1 = g(1) * t11 + g(2) * t9 + g(3) * t71;
t3 = g(1) * t12 + g(2) * t10 + g(3) * t69;
t50 = pkin(4) * t75 + t58 * t34 + t62;
t13 = t20 * t44 + t77;
t47 = pkin(4) * t73 + t61 + (-pkin(2) - t58) * t33;
t26 = qJ(6) * t69;
t25 = pkin(4) * t74;
t18 = t54 * t41;
t17 = t34 * t65 + t75;
t16 = -t63 + t74;
t15 = -t33 * t65 + t73;
t6 = t49 * t36;
t5 = t49 * t35;
t4 = g(1) * t10 - g(2) * t12;
t2 = [0, 0, 0, 0, 0, 0, t53, g(1) * t45 + g(2) * t42, 0, 0, 0, 0, 0, 0, 0, 0, t54, t20, 0, t53 * pkin(1), 0, 0, 0, 0, 0, 0, t54 * t44, -t18, -t20, -g(1) * (-t33 * pkin(2) + t61) - g(2) * t62, 0, 0, 0, 0, 0, 0, -g(1) * t15 - g(2) * t17, -g(1) * t14 - g(2) * t16, t18, -g(1) * t61 - g(2) * (t56 * t34 + t62) - (-pkin(2) - t56) * t80, 0, 0, 0, 0, 0, 0, t4, -t57, t18, -g(1) * t47 - g(2) * t50, 0, 0, 0, 0, 0, 0, t4, t18, t57, -g(1) * (-t10 * pkin(5) - t9 * qJ(6) + t47) - g(2) * (t12 * pkin(5) + t11 * qJ(6) + t50); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, t13, 0, 0, 0, 0, 0, 0, 0, 0, t49 * t43, -t49 * t40, -t13, -g(3) * t56 + t20 * (pkin(3) * t41 - pkin(8) * t44) 0, 0, 0, 0, 0, 0, t6, -t5, -t13, -g(3) * t58 + t20 * (t32 * t41 + t64) 0, 0, 0, 0, 0, 0, t6, -t13, t5, -g(3) * (t52 * t44 + t58) + t20 * (t64 - (-t32 - t52) * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t16 + g(2) * t14 + t40 * t77, g(1) * t17 - g(2) * t15 + t43 * t77, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, -g(1) * t25 + (g(2) * t72 + t13 * t40) * pkin(4), 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * (-pkin(4) * t63 + t25 + t59) - g(2) * (-t14 * pkin(4) + t60) - g(3) * (t26 + (-pkin(4) * t40 - pkin(5) * t35) * t41); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t3, 0, 0, 0, 0, 0, 0, 0, 0, t1, 0, -t3, -g(1) * t59 - g(2) * t60 - g(3) * (-pkin(5) * t71 + t26); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg  = t2;
