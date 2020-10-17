% Calculate inertial parameters regressor of gravitation load for
% S5RRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% g_base [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
% 
% Output:
% taug_reg [5x(5*10)]
%   inertial parameter regressor of gravitation joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:36
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function taug_reg = S5RRRPR11_gravloadJ_reg2_slag_vp(qJ, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR11_gravloadJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR11_gravloadJ_reg2_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR11_gravloadJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From gravload_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:35:14
% EndTime: 2019-12-31 21:35:15
% DurationCPUTime: 0.49s
% Computational Cost: add. (224->96), mult. (573->128), div. (0->0), fcn. (636->8), ass. (0->63)
t43 = cos(qJ(2));
t40 = sin(qJ(1));
t44 = cos(qJ(1));
t21 = g(1) * t44 + g(2) * t40;
t39 = sin(qJ(2));
t87 = t21 * t39;
t92 = g(3) * t43 - t87;
t42 = cos(qJ(3));
t38 = sin(qJ(3));
t68 = t44 * t38;
t18 = -t40 * t42 + t43 * t68;
t67 = t44 * t42;
t19 = t40 * t38 + t43 * t67;
t37 = sin(qJ(5));
t41 = cos(qJ(5));
t3 = t18 * t41 - t19 * t37;
t52 = t37 * t42 - t38 * t41;
t71 = t40 * t43;
t16 = t38 * t71 + t67;
t70 = t42 * t43;
t17 = t40 * t70 - t68;
t62 = -t16 * t41 + t17 * t37;
t89 = g(3) * t39;
t91 = -g(1) * t3 + g(2) * t62 + t52 * t89;
t4 = t18 * t37 + t19 * t41;
t51 = t37 * t38 + t41 * t42;
t53 = t16 * t37 + t17 * t41;
t86 = g(1) * t4 + g(2) * t53 + t51 * t89;
t83 = -pkin(3) - pkin(4);
t82 = g(1) * t40;
t32 = t39 * pkin(7);
t34 = t43 * pkin(2);
t75 = t38 * t39;
t74 = t39 * t40;
t73 = t39 * t42;
t72 = t39 * t44;
t69 = t43 * t44;
t66 = t34 + t32;
t65 = t44 * pkin(1) + t40 * pkin(6);
t64 = qJ(4) * t38;
t63 = -pkin(1) - t34;
t61 = -pkin(2) - t64;
t60 = -t16 * pkin(3) + t17 * qJ(4);
t59 = -t18 * pkin(3) + t19 * qJ(4);
t58 = pkin(3) * t70 + t43 * t64 + t66;
t57 = pkin(2) * t69 + pkin(7) * t72 + t65;
t35 = t44 * pkin(6);
t56 = -t17 * pkin(3) - t16 * qJ(4) + t35;
t55 = g(1) * t16 - g(2) * t18;
t54 = -g(2) * t44 + t82;
t48 = t19 * pkin(3) + t18 * qJ(4) + t57;
t47 = (t63 - t32) * t82;
t2 = g(1) * t18 + g(2) * t16 + g(3) * t75;
t46 = g(1) * t19 + g(2) * t17 + g(3) * t73;
t27 = pkin(7) * t69;
t24 = pkin(7) * t71;
t22 = qJ(4) * t73;
t20 = g(1) * t74 - g(2) * t72;
t8 = t21 * t43 + t89;
t7 = t92 * t42;
t6 = t92 * t38;
t5 = g(1) * t17 - g(2) * t19;
t1 = [0, 0, 0, 0, 0, 0, t54, t21, 0, 0, 0, 0, 0, 0, 0, 0, t54 * t43, -t20, -t21, -g(1) * (-t40 * pkin(1) + t35) - g(2) * t65, 0, 0, 0, 0, 0, 0, t5, -t55, t20, -g(1) * t35 - g(2) * t57 - t47, 0, 0, 0, 0, 0, 0, t5, t20, t55, -g(1) * t56 - g(2) * t48 - t47, 0, 0, 0, 0, 0, 0, g(1) * t53 - g(2) * t4, -g(1) * t62 - g(2) * t3, -t20, -g(1) * (-t17 * pkin(4) + t56) - g(2) * (t19 * pkin(4) - pkin(8) * t72 + t48) - ((-pkin(7) + pkin(8)) * t39 + t63) * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t92, t8, 0, 0, 0, 0, 0, 0, 0, 0, -t7, t6, -t8, -g(1) * (-pkin(2) * t72 + t27) - g(2) * (-pkin(2) * t74 + t24) - g(3) * t66, 0, 0, 0, 0, 0, 0, -t7, -t8, -t6, -g(1) * t27 - g(2) * t24 - g(3) * t58 + (pkin(3) * t42 - t61) * t87, 0, 0, 0, 0, 0, 0, -t92 * t51, t92 * t52, t8, -g(1) * (-pkin(8) * t69 + t27) - g(2) * (-pkin(8) * t71 + t24) - g(3) * (pkin(4) * t70 + t58) + (g(3) * pkin(8) + t21 * (-t83 * t42 - t61)) * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t2, t46, 0, 0, 0, 0, 0, 0, 0, 0, t2, 0, -t46, -g(1) * t59 - g(2) * t60 - g(3) * (-pkin(3) * t75 + t22), 0, 0, 0, 0, 0, 0, -t91, -t86, 0, -g(1) * (-t18 * pkin(4) + t59) - g(2) * (-t16 * pkin(4) + t60) - g(3) * (t83 * t75 + t22); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, t86, 0, 0;];
taug_reg = t1;
