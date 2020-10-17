% Calculate minimal parameter regressor of gravitation load for
% S6RRPRRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d5,theta3]';
% 
% Output:
% taug_reg [6x30]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 12:35
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6RRPRRP9_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP9_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RRPRRP9_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRPRRP9_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 18:34:30
% EndTime: 2019-05-06 18:34:31
% DurationCPUTime: 0.46s
% Computational Cost: add. (419->103), mult. (755->166), div. (0->0), fcn. (914->12), ass. (0->57)
t36 = sin(qJ(2));
t37 = sin(qJ(1));
t39 = cos(qJ(2));
t54 = cos(pkin(6));
t64 = cos(qJ(1));
t45 = t54 * t64;
t13 = t37 * t36 - t39 * t45;
t35 = sin(qJ(5));
t38 = cos(qJ(5));
t14 = t36 * t45 + t37 * t39;
t29 = pkin(11) + qJ(4);
t26 = sin(t29);
t27 = cos(t29);
t31 = sin(pkin(6));
t52 = t31 * t64;
t6 = t14 * t27 - t26 * t52;
t71 = -t13 * t38 + t6 * t35;
t70 = t13 * t35 + t6 * t38;
t59 = t31 * t36;
t12 = t54 * t26 + t27 * t59;
t49 = t37 * t54;
t16 = -t36 * t49 + t64 * t39;
t58 = t31 * t37;
t10 = t16 * t27 + t26 * t58;
t15 = t64 * t36 + t39 * t49;
t2 = -t10 * t35 + t15 * t38;
t56 = t38 * t39;
t69 = g(2) * t71 - g(3) * (-t12 * t35 - t31 * t56) - g(1) * t2;
t67 = g(2) * t13;
t66 = g(2) * t14;
t65 = g(3) * t31;
t61 = t27 * t35;
t60 = t27 * t38;
t57 = t35 * t39;
t55 = t64 * pkin(1) + pkin(8) * t58;
t30 = sin(pkin(11));
t53 = t30 * t58;
t51 = -t37 * pkin(1) + pkin(8) * t52;
t50 = pkin(5) * t35 + pkin(9) + qJ(3);
t48 = t30 * t52;
t5 = t14 * t26 + t27 * t52;
t9 = t16 * t26 - t27 * t58;
t47 = -g(1) * t5 + g(2) * t9;
t46 = g(1) * t13 - g(2) * t15;
t32 = cos(pkin(11));
t24 = t32 * pkin(3) + pkin(2);
t25 = t38 * pkin(5) + pkin(4);
t33 = -qJ(6) - pkin(10);
t44 = t25 * t27 - t26 * t33 + t24;
t11 = t26 * t59 - t54 * t27;
t43 = g(1) * t9 + g(2) * t5 + g(3) * t11;
t42 = g(1) * t10 + g(2) * t6 + g(3) * t12;
t4 = -g(1) * t15 + t39 * t65 - t67;
t41 = g(1) * t16 + g(3) * t59 + t66;
t3 = t10 * t38 + t15 * t35;
t1 = t4 * t26;
t7 = [0, g(1) * t37 - g(2) * t64, g(1) * t64 + g(2) * t37, 0, 0, 0, 0, 0, g(1) * t14 - g(2) * t16, -t46, -g(1) * (-t14 * t32 + t48) - g(2) * (t16 * t32 + t53) -g(1) * (t14 * t30 + t32 * t52) - g(2) * (-t16 * t30 + t32 * t58) t46, -g(1) * (-t14 * pkin(2) - t13 * qJ(3) + t51) - g(2) * (t16 * pkin(2) + t15 * qJ(3) + t55) 0, 0, 0, 0, 0, g(1) * t6 - g(2) * t10, t47, 0, 0, 0, 0, 0, g(1) * t70 - g(2) * t3, -g(1) * t71 - g(2) * t2, -t47, -g(1) * (pkin(3) * t48 - t50 * t13 - t14 * t24 - t25 * t6 + t33 * t5 + t51) - g(2) * (pkin(3) * t53 + t10 * t25 + t50 * t15 + t16 * t24 - t9 * t33 + t55); 0, 0, 0, 0, 0, 0, 0, 0, -t4, t41, -t4 * t32, t4 * t30, -t41, -g(1) * (-t15 * pkin(2) + t16 * qJ(3)) - g(2) * (-t13 * pkin(2) + t14 * qJ(3)) - (pkin(2) * t39 + qJ(3) * t36) * t65, 0, 0, 0, 0, 0, -t4 * t27, t1, 0, 0, 0, 0, 0, -g(1) * (-t15 * t60 + t16 * t35) - g(2) * (-t13 * t60 + t14 * t35) - (t27 * t56 + t35 * t36) * t65, -g(1) * (t15 * t61 + t16 * t38) - g(2) * (t13 * t61 + t14 * t38) - (-t27 * t57 + t36 * t38) * t65, -t1, -g(1) * (-t44 * t15 + t50 * t16) - t50 * t66 + t44 * t67 - (t50 * t36 + t44 * t39) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, t42, 0, 0, 0, 0, 0, t43 * t38, -t43 * t35, -t42, -g(1) * (-t10 * t33 - t9 * t25) - g(2) * (-t5 * t25 - t6 * t33) - g(3) * (-t11 * t25 - t12 * t33); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, g(1) * t3 + g(2) * t70 - g(3) * (-t12 * t38 + t31 * t57) 0, t69 * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43;];
taug_reg  = t7;
