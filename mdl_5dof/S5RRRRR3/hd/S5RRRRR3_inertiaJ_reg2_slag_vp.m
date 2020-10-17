% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRR3_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_inertiaJ_reg2_slag_vp: pkin has to be [5x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:56:29
% EndTime: 2019-12-05 18:56:33
% DurationCPUTime: 1.02s
% Computational Cost: add. (688->98), mult. (1591->210), div. (0->0), fcn. (1883->8), ass. (0->84)
t64 = sin(qJ(3));
t100 = t64 * pkin(1);
t54 = pkin(5) + t100;
t63 = sin(qJ(4));
t59 = t63 ^ 2;
t67 = cos(qJ(4));
t60 = t67 ^ 2;
t80 = t59 + t60;
t82 = t80 * t54;
t62 = sin(qJ(5));
t66 = cos(qJ(5));
t108 = -t62 * t63 + t66 * t67;
t65 = sin(qJ(2));
t68 = cos(qJ(3));
t69 = cos(qJ(2));
t41 = t64 * t65 - t68 * t69;
t36 = t41 ^ 2;
t107 = 0.2e1 * t41;
t44 = t64 * t69 + t68 * t65;
t106 = -0.2e1 * t44;
t56 = -t67 * pkin(3) - pkin(2);
t98 = t68 * pkin(1);
t46 = t56 - t98;
t105 = 0.2e1 * t46;
t104 = 0.2e1 * t56;
t102 = t41 * pkin(3);
t97 = t69 * pkin(1);
t18 = t41 * pkin(2) - t44 * pkin(5) - t97;
t15 = t67 * t18;
t7 = t15 + t102;
t93 = t63 * t18;
t4 = -t62 * t93 + t66 * t7;
t89 = t66 * t63;
t43 = t62 * t67 + t89;
t75 = t18 * t89;
t5 = t62 * t7 + t75;
t103 = t108 * t5 - t4 * t43;
t101 = t62 * pkin(3);
t99 = t66 * pkin(3);
t55 = -pkin(2) - t98;
t96 = pkin(2) - t55;
t38 = t44 ^ 2;
t95 = t59 * t38;
t31 = t63 * t41;
t92 = t63 * t44;
t91 = t63 * t54;
t90 = t63 * t67;
t32 = t67 * t41;
t87 = t67 * t44;
t86 = t67 * t54;
t27 = t43 * t54;
t28 = -t62 * t91 + t66 * t86;
t85 = t108 * t28 + t27 * t43;
t33 = t43 * pkin(5);
t34 = t108 * pkin(5);
t84 = t108 * t34 + t33 * t43;
t83 = t46 + t56;
t81 = t80 * pkin(5);
t79 = -0.2e1 * t97;
t78 = -0.2e1 * t31;
t77 = 0.2e1 * t32;
t76 = pkin(3) * t92;
t74 = -pkin(2) * t44 - pkin(5) * t41;
t73 = -t41 * t54 + t44 * t55;
t72 = pkin(1) ^ 2;
t71 = pkin(3) ^ 2;
t61 = t69 ^ 2;
t50 = 0.2e1 * t90;
t37 = t43 ^ 2;
t35 = t108 ^ 2;
t30 = t63 * t87;
t26 = t43 * t41;
t25 = t108 * t41;
t24 = 0.2e1 * t43 * t108;
t21 = t43 * t76;
t20 = t108 * t76;
t19 = (-t59 + t60) * t44;
t16 = (t108 * t62 - t43 * t66) * pkin(3);
t12 = -t62 * t92 + t66 * t87;
t10 = t43 * t44;
t9 = t12 * t43;
t8 = t10 * t108;
t3 = -t43 * t10 + t108 * t12;
t1 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t65 ^ 2, 0.2e1 * t65 * t69, 0, t61, 0, 0, 0, 0, 0, 0, t38, t41 * t106, 0, t36, 0, 0, t41 * t79, t44 * t79, 0, t61 * t72, t60 * t38, -0.2e1 * t38 * t90, t44 * t77, t95, t44 * t78, t36, t18 * t77, t18 * t78, t80 * t18 * t106, t80 * t18 ^ 2, t12 ^ 2, -0.2e1 * t12 * t10, t12 * t107, t10 ^ 2, -t10 * t107, t36, 0.2e1 * t10 * t76 + 0.2e1 * t4 * t41, 0.2e1 * t12 * t76 - 0.2e1 * t5 * t41, -0.2e1 * t5 * t10 - 0.2e1 * t4 * t12, t4 ^ 2 + t5 ^ 2 + t71 * t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, 0, t69, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t41, 0, 0, 0, (-t41 * t64 - t44 * t68) * pkin(1), 0, t30, t19, t31, -t30, t32, 0, t73 * t63, t73 * t67, 0, 0, t9, t3, t26, -t8, t25, 0, t46 * t10 - t27 * t41 - t20, t46 * t12 - t28 * t41 + t21, -t28 * t10 + t27 * t12 + t103, -t4 * t27 + t5 * t28 + t46 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t98, -0.2e1 * t100, 0, (t64 ^ 2 + t68 ^ 2) * t72, t59, t50, 0, t60, 0, 0, -0.2e1 * t55 * t67, 0.2e1 * t55 * t63, 0.2e1 * t82, t54 ^ 2 * t80 + t55 ^ 2, t37, t24, 0, t35, 0, 0, -t108 * t105, t43 * t105, 0.2e1 * t85, t27 ^ 2 + t28 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t44, 0, -t41, 0, 0, 0, 0, 0, t30, t19, t31, -t30, t32, 0, t74 * t63, t74 * t67, 0, 0, t9, t3, t26, -t8, t25, 0, t56 * t10 - t33 * t41 - t20, t56 * t12 - t34 * t41 + t21, -t34 * t10 + t33 * t12 + t103, -t4 * t33 + t5 * t34 + t56 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t98, -t100, 0, 0, t59, t50, 0, t60, 0, 0, t96 * t67, -t96 * t63, t81 + t82, -t55 * pkin(2) + pkin(5) * t82, t37, t24, 0, t35, 0, 0, -t83 * t108, t83 * t43, t84 + t85, t27 * t33 + t28 * t34 + t46 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t59, t50, 0, t60, 0, 0, 0.2e1 * pkin(2) * t67, -0.2e1 * pkin(2) * t63, 0.2e1 * t81, pkin(5) ^ 2 * t80 + pkin(2) ^ 2, t37, t24, 0, t35, 0, 0, -t108 * t104, t43 * t104, 0.2e1 * t84, t33 ^ 2 + t34 ^ 2 + t56 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t87, 0, -t92, t41, t15, -t93, 0, 0, 0, 0, t12, 0, -t10, t41, t41 * t99 + t4, -t75 + (-t7 - t102) * t62, (-t10 * t62 - t12 * t66) * pkin(3), (t4 * t66 + t5 * t62) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t67, 0, -t91, -t86, 0, 0, 0, 0, t43, 0, t108, 0, -t27, -t28, t16, (-t27 * t66 + t28 * t62) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t63, 0, t67, 0, -t63 * pkin(5), -t67 * pkin(5), 0, 0, 0, 0, t43, 0, t108, 0, -t33, -t34, t16, (-t33 * t66 + t34 * t62) * pkin(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t99, -0.2e1 * t101, 0, (t62 ^ 2 + t66 ^ 2) * t71; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t12, 0, -t10, t41, t4, -t5, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t108, 0, -t27, -t28, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, t108, 0, -t33, -t34, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t99, -t101, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0;];
MM_reg = t1;
