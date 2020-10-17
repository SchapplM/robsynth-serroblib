% Calculate inertial parameters regressor of joint inertia matrix for
% S5RRRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((5+1)*5/2)x(5*10)]
%   inertial parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S5RRRRP7_inertiaJ_reg2_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP7_inertiaJ_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP7_inertiaJ_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:57:44
% EndTime: 2019-12-31 21:57:48
% DurationCPUTime: 1.08s
% Computational Cost: add. (758->105), mult. (1503->188), div. (0->0), fcn. (1608->6), ass. (0->85)
t59 = sin(qJ(4));
t55 = t59 ^ 2;
t62 = cos(qJ(4));
t57 = t62 ^ 2;
t108 = t55 + t57;
t60 = sin(qJ(3));
t97 = t60 * pkin(2);
t48 = pkin(8) + t97;
t82 = t108 * t48;
t101 = -pkin(7) - pkin(6);
t64 = cos(qJ(2));
t39 = t101 * t64;
t63 = cos(qJ(3));
t61 = sin(qJ(2));
t75 = t101 * t61;
t16 = -t60 * t39 - t63 * t75;
t107 = t16 ^ 2;
t30 = t60 * t61 - t63 * t64;
t28 = t30 ^ 2;
t32 = t60 * t64 + t63 * t61;
t106 = 0.2e1 * t32;
t50 = -t64 * pkin(2) - pkin(1);
t105 = 0.2e1 * t50;
t104 = -0.2e1 * t59;
t103 = -0.2e1 * t62;
t102 = 0.2e1 * t64;
t100 = pkin(8) * t30;
t99 = t30 * pkin(4);
t98 = t59 * pkin(8);
t96 = t62 * pkin(8);
t95 = t63 * pkin(2);
t70 = pkin(4) * t59 - t62 * qJ(5);
t7 = t32 * t70 + t16;
t94 = t7 * t59;
t93 = t7 * t62;
t49 = -pkin(3) - t95;
t92 = pkin(3) - t49;
t10 = t30 * pkin(3) - t32 * pkin(8) + t50;
t18 = -t63 * t39 + t60 * t75;
t6 = t59 * t10 + t62 * t18;
t91 = t16 * t62;
t90 = t30 * t48;
t89 = t59 * t32;
t88 = t59 * t48;
t87 = t59 * t62;
t26 = t62 * t32;
t86 = t62 * t48;
t71 = -t62 * pkin(4) - t59 * qJ(5);
t35 = -pkin(3) + t71;
t27 = t35 - t95;
t85 = -t27 - t35;
t84 = t82 * pkin(8);
t83 = t108 * t48 ^ 2;
t81 = t108 * pkin(8) ^ 2;
t80 = t108 * pkin(8);
t56 = t61 ^ 2;
t58 = t64 ^ 2;
t79 = t56 + t58;
t78 = t30 * qJ(5);
t77 = t30 * t89;
t29 = t32 ^ 2;
t76 = t29 * t87;
t74 = -t62 * t10 + t59 * t18;
t73 = -pkin(3) * t32 - t100;
t3 = t78 + t6;
t4 = t74 - t99;
t1 = t3 * t62 + t4 * t59;
t2 = t59 * t74 + t6 * t62;
t72 = -t32 * t35 + t100;
t69 = t27 * t32 - t90;
t68 = t32 * t49 - t90;
t45 = -0.2e1 * t87;
t44 = 0.2e1 * t87;
t34 = 0.2e1 * t80;
t25 = t62 * t30;
t23 = t57 * t29;
t22 = t59 * t30;
t21 = t55 * t29;
t20 = t59 * t26;
t19 = 0.2e1 * t82;
t15 = t80 + t82;
t14 = t16 * t59;
t12 = 0.2e1 * t30 * t26;
t11 = (t55 - t57) * t32;
t5 = [0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t56, t61 * t102, 0, t58, 0, 0, pkin(1) * t102, -0.2e1 * pkin(1) * t61, 0.2e1 * t79 * pkin(6), t79 * pkin(6) ^ 2 + pkin(1) ^ 2, t29, -0.2e1 * t32 * t30, 0, t28, 0, 0, t30 * t105, t32 * t105, 0.2e1 * t16 * t32 - 0.2e1 * t18 * t30, t18 ^ 2 + t50 ^ 2 + t107, t23, -0.2e1 * t76, t12, t21, -0.2e1 * t77, t28, 0.2e1 * t16 * t89 - 0.2e1 * t30 * t74, 0.2e1 * t16 * t26 - 0.2e1 * t6 * t30, (-t59 * t6 + t62 * t74) * t106, t6 ^ 2 + t74 ^ 2 + t107, t23, t12, 0.2e1 * t76, t28, 0.2e1 * t77, t21, -0.2e1 * t4 * t30 + 0.2e1 * t7 * t89, (-t3 * t59 + t4 * t62) * t106, -0.2e1 * t26 * t7 + 0.2e1 * t3 * t30, t3 ^ 2 + t4 ^ 2 + t7 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t61, 0, t64, 0, -t61 * pkin(6), -t64 * pkin(6), 0, 0, 0, 0, t32, 0, -t30, 0, -t16, -t18, (-t30 * t60 - t32 * t63) * pkin(2), (-t16 * t63 + t18 * t60) * pkin(2), t20, -t11, t22, -t20, t25, 0, t59 * t68 - t91, t62 * t68 + t14, t2, t16 * t49 + t2 * t48, t20, t22, t11, 0, -t25, -t20, t59 * t69 - t93, t1, -t62 * t69 - t94, t1 * t48 + t7 * t27; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t95, -0.2e1 * t97, 0, (t60 ^ 2 + t63 ^ 2) * pkin(2) ^ 2, t55, t44, 0, t57, 0, 0, t49 * t103, 0.2e1 * t49 * t59, t19, t49 ^ 2 + t83, t55, 0, t45, 0, 0, t57, t27 * t103, t19, t27 * t104, t27 ^ 2 + t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, 0, -t30, 0, -t16, -t18, 0, 0, t20, -t11, t22, -t20, t25, 0, t59 * t73 - t91, t62 * t73 + t14, t2, -t16 * pkin(3) + pkin(8) * t2, t20, t22, t11, 0, -t25, -t20, -t59 * t72 - t93, t1, t62 * t72 - t94, pkin(8) * t1 + t7 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t95, -t97, 0, 0, t55, t44, 0, t57, 0, 0, t92 * t62, -t92 * t59, t15, -t49 * pkin(3) + t84, t55, 0, t45, 0, 0, t57, t85 * t62, t15, t85 * t59, t27 * t35 + t84; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, t55, t44, 0, t57, 0, 0, 0.2e1 * pkin(3) * t62, pkin(3) * t104, t34, pkin(3) ^ 2 + t81, t55, 0, t45, 0, 0, t57, t35 * t103, t34, t35 * t104, t35 ^ 2 + t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, 0, -t89, t30, -t74, -t6, 0, 0, 0, t26, 0, t30, t89, 0, -t74 + 0.2e1 * t99, t71 * t32, 0.2e1 * t78 + t6, -t4 * pkin(4) + t3 * qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t62, 0, -t88, -t86, 0, 0, 0, t59, 0, 0, -t62, 0, -t88, -t70, t86, -t70 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t62, 0, -t98, -t96, 0, 0, 0, t59, 0, 0, -t62, 0, -t98, -t70, t96, -t70 * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4), 0, 0.2e1 * qJ(5), pkin(4) ^ 2 + qJ(5) ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t30, t26, 0, t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t59, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -1, 0, 0, -pkin(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg = t5;
