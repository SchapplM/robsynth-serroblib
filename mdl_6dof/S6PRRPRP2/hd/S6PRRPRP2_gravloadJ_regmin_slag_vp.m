% Calculate minimal parameter regressor of gravitation load for
% S6PRRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1,theta4]';
% 
% Output:
% taug_reg [6x26]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 02:57
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S6PRRPRP2_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(3,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRPRP2_gravloadJ_regmin_slag_vp: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 02:54:48
% EndTime: 2021-01-16 02:54:51
% DurationCPUTime: 0.45s
% Computational Cost: add. (474->99), mult. (869->151), div. (0->0), fcn. (1067->12), ass. (0->67)
t55 = sin(qJ(3));
t58 = cos(qJ(3));
t82 = cos(pkin(6));
t52 = sin(pkin(6));
t56 = sin(qJ(2));
t86 = t52 * t56;
t101 = -t55 * t86 + t82 * t58;
t51 = sin(pkin(10));
t75 = t51 * t82;
t81 = cos(pkin(10));
t92 = cos(qJ(2));
t33 = t56 * t75 - t81 * t92;
t85 = t52 * t58;
t100 = t33 * t55 + t51 * t85;
t78 = t52 * t92;
t67 = t82 * t81;
t35 = t51 * t92 + t56 * t67;
t74 = t52 * t81;
t99 = -t35 * t55 - t58 * t74;
t50 = qJ(3) + pkin(11);
t48 = sin(t50);
t49 = cos(t50);
t98 = pkin(4) * t49 + pkin(9) * t48;
t87 = t51 * t52;
t62 = g(3) * (-t48 * t86 + t82 * t49) + g(2) * (-t35 * t48 - t49 * t74) + g(1) * (t33 * t48 + t49 * t87);
t54 = sin(qJ(5));
t89 = t49 * t54;
t57 = cos(qJ(5));
t88 = t49 * t57;
t36 = t81 * t56 + t92 * t75;
t47 = t58 * pkin(3) + pkin(2);
t53 = qJ(4) + pkin(8);
t84 = -t33 * t53 - t36 * t47;
t83 = t47 * t78 + t53 * t86;
t77 = t57 * t92;
t34 = t51 * t56 - t92 * t67;
t76 = -t34 * t47 + t35 * t53;
t15 = t33 * t49 - t48 * t87;
t72 = t54 * t78;
t71 = t100 * pkin(3);
t68 = t101 * pkin(3);
t27 = t82 * t48 + t49 * t86;
t20 = t27 * t54 + t52 * t77;
t17 = t35 * t49 - t48 * t74;
t6 = t17 * t54 - t34 * t57;
t8 = -t15 * t54 - t36 * t57;
t1 = g(1) * t8 + g(2) * t6 + g(3) * t20;
t21 = t27 * t57 - t72;
t7 = t17 * t57 + t34 * t54;
t9 = -t15 * t57 + t36 * t54;
t65 = g(1) * t9 + g(2) * t7 + g(3) * t21;
t10 = -t34 * t89 - t35 * t57;
t12 = t33 * t57 - t36 * t89;
t22 = t49 * t72 - t57 * t86;
t64 = g(1) * t12 + g(2) * t10 + g(3) * t22;
t63 = g(1) * t15 - g(2) * t17 - g(3) * t27;
t61 = -g(1) * t33 + g(2) * t35 + g(3) * t86;
t60 = g(2) * t99;
t59 = g(1) * t36 + g(2) * t34 - g(3) * t78;
t23 = (t49 * t77 + t54 * t56) * t52;
t13 = -t33 * t54 - t36 * t88;
t11 = -t34 * t88 + t35 * t54;
t5 = t59 * t48;
t4 = t62 * t57;
t3 = t62 * t54;
t2 = -g(1) * t13 - g(2) * t11 - g(3) * t23;
t14 = [-g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(3); 0, 0, t59, t61, 0, 0, 0, 0, 0, t59 * t58, -t59 * t55, t59 * t49, -t5, -t61, -g(1) * t84 - g(2) * t76 - g(3) * t83, 0, 0, 0, 0, 0, t2, t64, t2, t5, -t64, -g(1) * (t13 * pkin(5) + t12 * qJ(6) - t36 * t98 + t84) - g(2) * (t11 * pkin(5) + t10 * qJ(6) - t34 * t98 + t76) - g(3) * (t23 * pkin(5) + t22 * qJ(6) + t98 * t78 + t83); 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * t100 - g(3) * t101 - t60, -g(1) * (t33 * t58 - t55 * t87) - g(2) * (-t35 * t58 + t55 * t74) - g(3) * (-t82 * t55 - t56 * t85), -t62, -t63, 0, -pkin(3) * t60 - g(1) * t71 - g(3) * t68, 0, 0, 0, 0, 0, -t4, t3, -t4, t63, -t3, -g(1) * (-t15 * pkin(9) + t71) - g(2) * (pkin(3) * t99 + t17 * pkin(9)) - g(3) * (t27 * pkin(9) + t68) - t62 * (pkin(5) * t57 + qJ(6) * t54 + pkin(4)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t59; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1, t65, t1, 0, -t65, -g(1) * (-t8 * pkin(5) + t9 * qJ(6)) - g(2) * (-t6 * pkin(5) + t7 * qJ(6)) - g(3) * (-t20 * pkin(5) + t21 * qJ(6)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t1;];
taug_reg = t14;
