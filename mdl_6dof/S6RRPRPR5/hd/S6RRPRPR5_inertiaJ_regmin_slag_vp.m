% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d4,d6,theta3,theta5]';
% 
% Output:
% MM_reg [((6+1)*6/2)x30]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:37
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRPRPR5_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPR5_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRPRPR5_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 14:02:25
% EndTime: 2019-05-06 14:02:30
% DurationCPUTime: 1.18s
% Computational Cost: add. (2141->165), mult. (5241->327), div. (0->0), fcn. (6033->12), ass. (0->106)
t75 = sin(pkin(11));
t76 = sin(pkin(6));
t78 = cos(pkin(11));
t82 = sin(qJ(2));
t85 = cos(qJ(2));
t46 = (t75 * t85 + t78 * t82) * t76;
t79 = cos(pkin(6));
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t35 = t46 * t81 - t79 * t84;
t123 = -0.2e1 * t35;
t77 = cos(pkin(12));
t68 = -t77 * pkin(5) - pkin(4);
t122 = 0.2e1 * t68;
t121 = 0.2e1 * t81;
t120 = -0.2e1 * t84;
t119 = 0.2e1 * t84;
t118 = pkin(1) * t82;
t117 = t84 * pkin(4);
t108 = t76 * t82;
t63 = t79 * t85 * pkin(1);
t98 = pkin(8) + qJ(3);
t39 = t79 * pkin(2) - t98 * t108 + t63;
t107 = t76 * t85;
t93 = t79 * t118;
t43 = t98 * t107 + t93;
t24 = t78 * t39 - t75 * t43;
t21 = -t79 * pkin(3) - t24;
t36 = t46 * t84 + t79 * t81;
t13 = t35 * pkin(4) - t36 * qJ(5) + t21;
t74 = sin(pkin(12));
t25 = t75 * t39 + t78 * t43;
t22 = t79 * pkin(9) + t25;
t45 = -t78 * t107 + t75 * t108;
t57 = (-pkin(2) * t85 - pkin(1)) * t76;
t28 = t45 * pkin(3) - t46 * pkin(9) + t57;
t15 = t84 * t22 + t81 * t28;
t9 = t45 * qJ(5) + t15;
t6 = t74 * t13 + t77 * t9;
t14 = -t81 * t22 + t84 * t28;
t10 = -t45 * pkin(4) - t14;
t116 = t10 * t74;
t115 = t10 * t77;
t26 = t36 * t74 - t45 * t77;
t27 = t36 * t77 + t45 * t74;
t80 = sin(qJ(6));
t83 = cos(qJ(6));
t17 = -t80 * t26 + t83 * t27;
t114 = t17 * t84;
t55 = t80 * t74 - t83 * t77;
t48 = t55 * t81;
t113 = t48 * t35;
t56 = t83 * t74 + t80 * t77;
t112 = t56 * t84;
t66 = t75 * pkin(2) + pkin(9);
t111 = t66 * t74;
t70 = t76 ^ 2;
t110 = t70 * t85;
t109 = t74 * t81;
t106 = t77 * t81;
t105 = t81 * t35;
t104 = t81 * t45;
t61 = t81 * t66;
t16 = t83 * t26 + t80 * t27;
t103 = t84 * t16;
t102 = t84 * t55;
t101 = t84 * t66;
t100 = t84 * t74;
t99 = t84 * t77;
t97 = pkin(10) + qJ(5);
t67 = -t78 * pkin(2) - pkin(3);
t54 = -t81 * qJ(5) - t117 + t67;
t33 = t74 * t54 + t66 * t99;
t96 = t74 ^ 2 + t77 ^ 2;
t95 = qJ(5) * t35;
t94 = 0.2e1 * t76 * t79;
t5 = t77 * t13 - t74 * t9;
t92 = t96 * qJ(5);
t91 = -t5 * t74 + t6 * t77;
t90 = -pkin(4) * t81 + qJ(5) * t84;
t89 = -t26 * t77 + t27 * t74;
t50 = t77 * t54;
t32 = -t66 * t100 + t50;
t88 = -t32 * t74 + t33 * t77;
t73 = t84 ^ 2;
t72 = t81 ^ 2;
t60 = t97 * t77;
t59 = t97 * t74;
t53 = pkin(8) * t107 + t93;
t52 = -pkin(8) * t108 + t63;
t51 = pkin(5) * t109 + t61;
t47 = t56 * t81;
t44 = t84 * t45;
t38 = -t80 * t59 + t83 * t60;
t37 = -t83 * t59 - t80 * t60;
t31 = -pkin(10) * t109 + t33;
t30 = -pkin(10) * t106 + t50 + (-pkin(5) - t111) * t84;
t29 = t47 * t35;
t19 = t80 * t30 + t83 * t31;
t18 = t83 * t30 - t80 * t31;
t7 = t26 * pkin(5) + t10;
t4 = -t26 * pkin(10) + t6;
t3 = t35 * pkin(5) - t27 * pkin(10) + t5;
t2 = t80 * t3 + t83 * t4;
t1 = t83 * t3 - t80 * t4;
t8 = [1, 0, 0, t70 * t82 ^ 2, 0.2e1 * t82 * t110, t82 * t94, t85 * t94, t79 ^ 2, 0.2e1 * pkin(1) * t110 + 0.2e1 * t52 * t79, -0.2e1 * t70 * t118 - 0.2e1 * t53 * t79, -0.2e1 * t24 * t46 - 0.2e1 * t25 * t45, t24 ^ 2 + t25 ^ 2 + t57 ^ 2, t36 ^ 2, t36 * t123, 0.2e1 * t36 * t45, t45 * t123, t45 ^ 2, 0.2e1 * t14 * t45 + 0.2e1 * t21 * t35, -0.2e1 * t15 * t45 + 0.2e1 * t21 * t36, 0.2e1 * t10 * t26 + 0.2e1 * t5 * t35, 0.2e1 * t10 * t27 - 0.2e1 * t6 * t35, -0.2e1 * t6 * t26 - 0.2e1 * t5 * t27, t10 ^ 2 + t5 ^ 2 + t6 ^ 2, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t35, t16 * t123, t35 ^ 2, 0.2e1 * t1 * t35 + 0.2e1 * t7 * t16, 0.2e1 * t7 * t17 - 0.2e1 * t2 * t35; 0, 0, 0, 0, 0, t108, t107, t79, t52, -t53 (-t45 * t75 - t46 * t78) * pkin(2) (t24 * t78 + t25 * t75) * pkin(2), t36 * t81, t36 * t84 - t105, t104, t44, 0, -t21 * t84 + t67 * t35 - t45 * t61, -t45 * t101 + t21 * t81 + t67 * t36, t32 * t35 - t5 * t84 + (t26 * t66 + t116) * t81, -t33 * t35 + t6 * t84 + (t27 * t66 + t115) * t81, -t33 * t26 - t32 * t27 + (-t5 * t77 - t6 * t74) * t81, t10 * t61 + t5 * t32 + t6 * t33, -t17 * t48, t48 * t16 - t17 * t47, -t113 - t114, -t29 + t103, -t35 * t84, -t1 * t84 + t51 * t16 + t18 * t35 + t7 * t47, t51 * t17 - t19 * t35 + t2 * t84 - t7 * t48; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0 (t75 ^ 2 + t78 ^ 2) * pkin(2) ^ 2, t72, t81 * t119, 0, 0, 0, t67 * t120, t67 * t121, 0.2e1 * t72 * t111 - 0.2e1 * t32 * t84, 0.2e1 * t72 * t66 * t77 + 0.2e1 * t33 * t84 (-t32 * t77 - t33 * t74) * t121, t72 * t66 ^ 2 + t32 ^ 2 + t33 ^ 2, t48 ^ 2, 0.2e1 * t48 * t47, -t48 * t120, t47 * t119, t73, -0.2e1 * t18 * t84 + 0.2e1 * t51 * t47, 0.2e1 * t19 * t84 - 0.2e1 * t51 * t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, 0, 0, 0, 0, 0, t44, -t104, -t74 * t105 - t84 * t26, -t77 * t105 - t84 * t27, t89 * t81, -t10 * t84 + t81 * t91, 0, 0, 0, 0, 0, -t29 - t103, t113 - t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0 (t88 - t101) * t81, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t72 * t96 + t73, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t36, -t35, t45, t14, -t15, -pkin(4) * t26 - t74 * t95 - t115, -pkin(4) * t27 - t77 * t95 + t116, t89 * qJ(5) + t91, -t10 * pkin(4) + qJ(5) * t91, t17 * t56, -t56 * t16 - t17 * t55, t56 * t35, -t55 * t35, 0, t68 * t16 + t37 * t35 + t7 * t55, t68 * t17 - t38 * t35 + t7 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t84, 0, -t61, -t101, -t77 * t61 + t90 * t74, t74 * t61 + t90 * t77, t88, -pkin(4) * t61 + qJ(5) * t88, -t48 * t56, -t56 * t47 + t48 * t55, -t112, t102, 0, -t37 * t84 + t68 * t47 + t51 * t55, t38 * t84 - t48 * t68 + t51 * t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t81, t99, -t100, t96 * t81, t81 * t92 + t117, 0, 0, 0, 0, 0, -t102, -t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0.2e1 * pkin(4) * t77, -0.2e1 * pkin(4) * t74, 0.2e1 * t92, t96 * qJ(5) ^ 2 + pkin(4) ^ 2, t56 ^ 2, -0.2e1 * t56 * t55, 0, 0, 0, t55 * t122, t56 * t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t26, t27, 0, t10, 0, 0, 0, 0, 0, t16, t17; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t109, t106, 0, t61, 0, 0, 0, 0, 0, t47, -t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t77, t74, 0, -pkin(4), 0, 0, 0, 0, 0, t55, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t35, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, -t47, -t84, t18, -t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t47, t48; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t56, -t55, 0, t37, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t8;
