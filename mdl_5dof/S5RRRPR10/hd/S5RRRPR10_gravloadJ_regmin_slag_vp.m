% Calculate minimal parameter regressor of gravitation load for
% S5RRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5,theta4]';
% 
% Output:
% taug_reg [5x28]
%   minimal parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 23:43
% Revision: 24b2e7d74a0c1a3b64fa2f8f5ad758691ad61af3 (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR10_gravloadJ_regmin_slag_vp(qJ, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(10,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5RRRPR10_gravloadJ_regmin_slag_vp: pkin has to be [10x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 23:41:17
% EndTime: 2021-01-15 23:41:19
% DurationCPUTime: 0.50s
% Computational Cost: add. (276->88), mult. (545->154), div. (0->0), fcn. (655->12), ass. (0->63)
t30 = sin(qJ(2));
t31 = sin(qJ(1));
t57 = cos(pkin(5));
t65 = cos(qJ(2));
t47 = t57 * t65;
t66 = cos(qJ(1));
t14 = t31 * t30 - t66 * t47;
t34 = -t66 * t30 - t31 * t47;
t26 = sin(pkin(5));
t68 = g(3) * t26;
t36 = -g(1) * t34 + g(2) * t14 - t65 * t68;
t49 = t66 * t65;
t53 = t30 * t57;
t12 = t31 * t53 - t49;
t25 = qJ(3) + pkin(10);
t23 = sin(t25);
t24 = cos(t25);
t61 = t26 * t31;
t2 = t12 * t23 + t24 * t61;
t54 = t31 * t65;
t13 = -t66 * t53 - t54;
t56 = t26 * t66;
t40 = t13 * t23 - t24 * t56;
t62 = t26 * t30;
t73 = g(2) * t40 + g(3) * (-t23 * t62 + t57 * t24) + g(1) * t2;
t17 = t23 * t56;
t41 = t13 * t24 + t17;
t44 = -t12 * t24 + t23 * t61;
t7 = t57 * t23 + t24 * t62;
t75 = g(1) * t44 - g(2) * t41 + g(3) * t7;
t32 = cos(qJ(5));
t28 = sin(qJ(5));
t64 = t34 * t28;
t74 = t44 * t32 - t64;
t72 = g(1) * t41 + g(2) * t44;
t29 = sin(qJ(3));
t67 = t29 * pkin(3);
t63 = t24 * t32;
t8 = t28 * t14;
t60 = t29 * t26;
t59 = t32 * t13;
t33 = cos(qJ(3));
t58 = t33 * t26;
t55 = t28 * t65;
t52 = t57 * t29;
t51 = t57 * t33;
t48 = g(1) * t14 + g(2) * t34;
t46 = t32 * t17 + t24 * t59 - t8;
t45 = t29 * t49;
t43 = -t12 * t33 + t31 * t60;
t22 = t33 * pkin(3) + pkin(2);
t27 = qJ(4) + pkin(8);
t42 = t22 * t65 + t27 * t30;
t39 = t30 * t52 + t58;
t38 = -g(1) * t12 - g(2) * t13 + g(3) * t62;
t37 = g(3) * (-t30 * t60 + t51);
t11 = -t22 * t30 + t27 * t65;
t10 = pkin(1) + t42;
t9 = t39 * pkin(3);
t5 = t22 * t47 + t27 * t53;
t4 = -t13 * t29 + t33 * t56;
t3 = t27 * t47 - t22 * t53 + t26 * (pkin(7) + t67);
t1 = [0, g(1) * t31 - g(2) * t66, g(1) * t66 + g(2) * t31, 0, 0, 0, 0, 0, -g(1) * t13 + g(2) * t12, -t48, 0, 0, 0, 0, 0, -g(1) * ((-t30 * t51 + t60) * t66 - t33 * t54) - g(2) * t43, -g(1) * t4 - g(2) * (t12 * t29 + t31 * t58), -t72, g(1) * t40 - g(2) * t2, t48, -g(1) * (-t10 * t31 + t3 * t66) - g(2) * (t10 * t66 + t3 * t31), 0, 0, 0, 0, 0, -g(1) * t46 - g(2) * t74, t72 * t28 + t48 * t32; 0, 0, 0, 0, 0, 0, 0, 0, t36, t38, 0, 0, 0, 0, 0, t36 * t33, -t36 * t29, t36 * t24, -t36 * t23, -t38, -g(1) * (t11 * t66 - t5 * t31) - g(2) * (t11 * t31 + t5 * t66) - t42 * t68, 0, 0, 0, 0, 0, -g(1) * (-t12 * t28 + t34 * t63) - g(2) * (-t28 * t13 - t14 * t63) - (t28 * t30 + t65 * t63) * t68, -g(1) * (-t12 * t32 - t24 * t64) - g(2) * (t24 * t8 - t59) - (-t24 * t55 + t30 * t32) * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -g(1) * (t39 * t31 - t45) + g(2) * t4 - t37, g(1) * t43 - g(2) * (t13 * t33 + t29 * t56) - g(3) * (-t30 * t58 - t52), -t73, t75, 0, -g(1) * (-pkin(3) * t45 + t9 * t31) - g(2) * (-t54 * t67 - t9 * t66) - pkin(3) * t37, 0, 0, 0, 0, 0, -t73 * t32, t73 * t28; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t36, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75 * t28 - t36 * t32, g(1) * t74 - g(2) * t46 - g(3) * (t26 * t55 - t7 * t32);];
taug_reg = t1;
