% Calculate minimal parameter regressor of gravitation load for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-16 00:37
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRRP10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_gravloadJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-16 00:35:52
% EndTime: 2021-01-16 00:35:55
% DurationCPUTime: 0.58s
% Computational Cost: add. (335->97), mult. (791->173), div. (0->0), fcn. (946->10), ass. (0->75)
t47 = cos(pkin(5));
t46 = sin(pkin(5));
t54 = cos(qJ(3));
t74 = t54 * t46;
t50 = sin(qJ(3));
t51 = sin(qJ(2));
t80 = t50 * t51;
t25 = t47 * t80 + t74;
t56 = cos(qJ(1));
t52 = sin(qJ(1));
t55 = cos(qJ(2));
t70 = t56 * t55;
t69 = g(1) * (-t25 * t52 + t50 * t70) + g(3) * (t46 * t80 - t47 * t54);
t77 = t52 * t55;
t94 = t69 + g(2) * (t25 * t56 + t50 * t77);
t71 = t56 * t51;
t33 = t47 * t77 + t71;
t78 = t52 * t51;
t32 = -t47 * t70 + t78;
t85 = g(3) * t46;
t65 = g(2) * t32 - t55 * t85;
t58 = -g(1) * t33 - t65;
t92 = t58 * t50;
t31 = -t47 * t71 - t77;
t19 = -t31 * t50 + t56 * t74;
t30 = t47 * t78 - t70;
t91 = g(2) * (t30 * t50 + t52 * t74) + g(1) * t19;
t81 = t50 * t46;
t20 = -t30 * t54 + t52 * t81;
t83 = t47 * t50;
t27 = t51 * t74 + t83;
t90 = g(2) * (t31 * t54 + t56 * t81) - g(3) * t27 - g(1) * t20;
t84 = t46 * t51;
t49 = sin(qJ(4));
t82 = t49 * t55;
t79 = t51 * t54;
t53 = cos(qJ(4));
t76 = t53 * t51;
t75 = t53 * t55;
t73 = t54 * t55;
t36 = t49 * t51 + t53 * t73;
t72 = t56 * t36;
t45 = pkin(4) * t53 + pkin(3);
t48 = qJ(5) + pkin(9);
t29 = -t45 * t50 + t48 * t54;
t28 = t47 * t79 - t81;
t16 = t28 * t49 + t47 * t75;
t35 = t49 * t73 - t76;
t68 = t16 * t56 + t35 * t52;
t66 = t45 * t54 + t48 * t50;
t24 = pkin(2) + t66;
t44 = pkin(4) * t49 + pkin(8);
t67 = t24 * t55 + t44 * t51;
t15 = -t24 * t51 + t44 * t55;
t62 = -t28 * t53 + t47 * t82;
t59 = t49 * t79 + t75;
t23 = t35 * pkin(4);
t22 = t33 * t54;
t21 = t35 * t47;
t18 = -t28 * t56 - t52 * t73;
t14 = pkin(1) + t67;
t13 = (t47 * t59 - t49 * t81) * pkin(4);
t12 = t67 * t47;
t11 = t29 * t47 * t51 - t46 * t66;
t10 = -t52 * t62 - t72;
t9 = -t46 * (-pkin(7) + t29) + t15 * t47;
t8 = t94 * t53;
t7 = t94 * t49;
t6 = -g(1) * (t18 * t53 - t32 * t49) + g(2) * t10;
t5 = -g(1) * t68 - g(2) * (t16 * t52 - t56 * t35);
t4 = -g(1) * (-t22 * t53 - t30 * t49) - g(2) * (t52 * (-t54 * t76 + t82) + t47 * t72) - t36 * t85;
t3 = -g(1) * (t21 * t52 + t56 * t59) - g(2) * (-t21 * t56 + t52 * t59) + t35 * t85;
t2 = -g(1) * (-(-t28 * t52 + t54 * t70) * t49 + t33 * t53) + g(2) * t68 - g(3) * (-t27 * t49 - t46 * t75);
t1 = -g(1) * t10 - g(2) * (-t52 * t36 + t56 * t62) - g(3) * (-t27 * t53 + t46 * t82);
t17 = [0, g(1) * t52 - g(2) * t56, g(1) * t56 + g(2) * t52, 0, 0, 0, 0, 0, -g(1) * t31 + g(2) * t30, -g(1) * t32 + g(2) * t33, 0, 0, 0, 0, 0, -g(1) * t18 - g(2) * t20, -t91, 0, 0, 0, 0, 0, t6, t5, t6, t5, t91, -g(1) * (-t14 * t52 + t56 * t9) - g(2) * (t14 * t56 + t52 * t9); 0, 0, 0, 0, 0, 0, 0, 0, -t58, -g(1) * t30 - g(2) * t31 + g(3) * t84, 0, 0, 0, 0, 0, g(1) * t22 + t54 * t65, t92, 0, 0, 0, 0, 0, t4, t3, t4, t3, -t92, -g(1) * (-t12 * t52 + t15 * t56) - g(2) * (t12 * t56 + t15 * t52) - t67 * t85; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, g(2) * t19 + t69, -t90, 0, 0, 0, 0, 0, t8, -t7, t8, -t7, t90, -g(1) * (-t11 * t52 + t29 * t70) - g(2) * (t11 * t56 + t29 * t77) - g(3) * (t29 * t84 + t47 * t66); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t1, t2, t1, 0, -g(1) * (t13 * t52 - t23 * t56) - g(2) * (-t13 * t56 - t23 * t52) - g(3) * (-t46 * t59 - t49 * t83) * pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t94;];
taug_reg = t17;
