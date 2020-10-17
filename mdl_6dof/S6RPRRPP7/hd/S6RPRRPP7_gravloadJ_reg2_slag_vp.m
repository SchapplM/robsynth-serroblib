% Calculate inertial parameters regressor of gravitation load for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
% 
% Output:
% taug_reg [6x(6*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RPRRPP7_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_gravloadJ_reg2_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRPP7_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 21:48:31
% EndTime: 2019-05-05 21:48:32
% DurationCPUTime: 0.38s
% Computational Cost: add. (202->97), mult. (492->116), div. (0->0), fcn. (515->6), ass. (0->60)
t39 = sin(qJ(3));
t38 = sin(qJ(4));
t43 = cos(qJ(1));
t66 = t43 * t38;
t40 = sin(qJ(1));
t41 = cos(qJ(4));
t69 = t40 * t41;
t16 = t39 * t66 + t69;
t65 = t43 * t41;
t70 = t40 * t38;
t17 = t39 * t65 - t70;
t81 = t17 * pkin(4) + t16 * qJ(5);
t80 = -pkin(1) - pkin(7);
t79 = -pkin(4) - pkin(5);
t78 = pkin(8) * t39;
t77 = g(1) * t40;
t42 = cos(qJ(3));
t76 = g(2) * t42;
t75 = g(2) * t43;
t74 = g(3) * t39;
t34 = t42 * pkin(8);
t73 = t38 * t42;
t72 = t39 * t40;
t71 = t39 * t43;
t68 = t40 * t42;
t67 = t41 * t42;
t64 = -pkin(8) + qJ(6);
t63 = pkin(3) * t68 + pkin(8) * t72;
t62 = t43 * pkin(1) + t40 * qJ(2);
t60 = t42 * qJ(6);
t59 = pkin(8) * t68;
t58 = t38 * t68;
t57 = t40 * t67;
t56 = t43 * pkin(7) + t62;
t55 = t80 * t40;
t54 = -qJ(5) * t38 - pkin(3);
t14 = t39 * t70 - t65;
t15 = t39 * t69 + t66;
t53 = -t14 * pkin(4) + t15 * qJ(5);
t52 = t16 * pkin(4) - t17 * qJ(5);
t51 = pkin(4) * t57 + qJ(5) * t58 + t63;
t50 = pkin(3) * t72 + t56;
t49 = g(1) * t16 + g(2) * t14;
t23 = g(1) * t43 + g(2) * t40;
t22 = -t75 + t77;
t33 = t43 * qJ(2);
t48 = pkin(3) * t71 - t43 * t34 + t33;
t47 = -pkin(4) * t41 + t54;
t46 = t15 * pkin(4) + t14 * qJ(5) + t50;
t2 = g(1) * t14 - g(2) * t16 + g(3) * t73;
t45 = g(1) * t15 - g(2) * t17 + g(3) * t67;
t44 = t55 + t48;
t9 = -t22 * t42 + t74;
t24 = qJ(5) * t67;
t18 = t23 * t42;
t8 = g(1) * t72 - g(2) * t71 + g(3) * t42;
t7 = t9 * t41;
t6 = g(1) * t58 - t38 * t74 - t66 * t76;
t5 = -g(1) * t17 - g(2) * t15;
t1 = [0, 0, 0, 0, 0, 0, t22, t23, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, -t23, -g(1) * (-t40 * pkin(1) + t33) - g(2) * t62, 0, 0, 0, 0, 0, 0, -t23 * t39, -t18, t22, -g(1) * (t33 + t55) - g(2) * t56, 0, 0, 0, 0, 0, 0, t5, t49, t18, -g(1) * t44 - g(2) * (t50 - t59) 0, 0, 0, 0, 0, 0, t5, t18, -t49, -g(1) * (t44 + t81) - g(2) * (t46 - t59) 0, 0, 0, 0, 0, 0, t5, -t49, -t18, -g(1) * (t17 * pkin(5) + t43 * t60 + t48 + t81) - g(2) * (t15 * pkin(5) + t46) + (-g(1) * t80 - t64 * t76) * t40; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t22; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9, t8, 0, 0, 0, 0, 0, 0, 0, 0, t7, t6, -t8, -g(1) * t63 - g(3) * (-t39 * pkin(3) + t34) - (-pkin(3) * t42 - t78) * t75, 0, 0, 0, 0, 0, 0, t7, -t8, -t6, -g(1) * t51 - g(3) * t34 - t47 * t74 - (t47 * t42 - t78) * t75, 0, 0, 0, 0, 0, 0, t7, -t6, t8, -g(1) * (pkin(5) * t57 + t51) - g(3) * (t34 - t60) + (qJ(6) * t77 - g(3) * (-pkin(5) * t41 + t47)) * t39 - (t64 * t39 + (t79 * t41 + t54) * t42) * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t45, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t45, -g(1) * t53 - g(2) * t52 - g(3) * (-pkin(4) * t73 + t24) 0, 0, 0, 0, 0, 0, t2, -t45, 0, -g(1) * (-t14 * pkin(5) + t53) - g(2) * (t16 * pkin(5) + t52) - g(3) * (t79 * t73 + t24); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t9;];
taug_reg  = t1;
