% Calculate minimal parameter regressor of joint inertia matrix for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RRRRPP3_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_inertiaJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-07 18:10:48
% EndTime: 2019-05-07 18:10:51
% DurationCPUTime: 1.12s
% Computational Cost: add. (1099->132), mult. (1982->227), div. (0->0), fcn. (2181->6), ass. (0->97)
t69 = sin(qJ(3));
t109 = t69 * pkin(2);
t55 = pkin(9) + t109;
t68 = sin(qJ(4));
t64 = t68 ^ 2;
t71 = cos(qJ(4));
t65 = t71 ^ 2;
t90 = t64 + t65;
t92 = t90 * t55;
t58 = qJ(5) * t71;
t122 = -pkin(4) * t68 + t58;
t66 = pkin(4) + qJ(6);
t89 = t68 * qJ(5);
t121 = -t66 * t71 - t89;
t104 = cos(qJ(3));
t105 = cos(qJ(2));
t70 = sin(qJ(2));
t37 = t104 * t70 + t69 * t105;
t120 = 0.2e1 * t37;
t57 = -t105 * pkin(2) - pkin(1);
t119 = 0.2e1 * t57;
t118 = 0.2e1 * t66;
t117 = -0.2e1 * t68;
t116 = -0.2e1 * t71;
t115 = 0.2e1 * t71;
t36 = -t104 * t105 + t69 * t70;
t29 = t71 * t37;
t15 = t36 * pkin(3) - t37 * pkin(9) + t57;
t45 = (-pkin(8) - pkin(7)) * t70;
t85 = t105 * pkin(7);
t47 = t105 * pkin(8) + t85;
t21 = t104 * t47 + t69 * t45;
t97 = -t71 * t15 + t68 * t21;
t75 = -pkin(5) * t29 - t97;
t3 = -t66 * t36 - t75;
t32 = t36 * qJ(5);
t10 = t68 * t15 + t71 * t21;
t99 = t68 * t37;
t76 = -pkin(5) * t99 + t10;
t5 = t32 + t76;
t114 = t3 * t68 + t5 * t71;
t112 = pkin(9) * t36;
t111 = t36 * pkin(4);
t110 = t68 * pkin(5);
t20 = -t104 * t45 + t69 * t47;
t83 = pkin(4) * t99 + t20;
t8 = (qJ(6) * t68 - t58) * t37 + t83;
t108 = t8 * t68;
t107 = t8 * t71;
t84 = t104 * pkin(2);
t56 = -t84 - pkin(3);
t106 = pkin(3) - t56;
t12 = -t37 * t58 + t83;
t103 = t12 * t68;
t102 = t20 * t71;
t101 = t36 * t55;
t98 = t68 * t71;
t30 = -pkin(3) + t121;
t22 = -t84 + t30;
t96 = -t22 - t30;
t49 = t68 * t55;
t33 = t49 + t110;
t51 = t71 * t55;
t63 = t71 * pkin(5);
t34 = t51 + t63;
t95 = t33 * t68 + t34 * t71;
t80 = -t71 * pkin(4) - t89;
t42 = -pkin(3) + t80;
t31 = -t84 + t42;
t94 = t31 + t42;
t60 = t68 * pkin(9);
t44 = t60 + t110;
t62 = t71 * pkin(9);
t46 = t62 + t63;
t93 = t44 * t68 + t46 * t71;
t91 = t90 * pkin(9);
t88 = 0.2e1 * t105;
t87 = -0.2e1 * t37 * t36;
t86 = t32 + t10;
t82 = -pkin(3) * t37 - t112;
t7 = t97 - t111;
t1 = t7 * t68 + t71 * t86;
t81 = -t37 * t42 + t112;
t78 = -t31 * t37 + t101;
t77 = t37 * t56 - t101;
t73 = qJ(5) ^ 2;
t72 = 0.2e1 * qJ(5);
t52 = 0.2e1 * t98;
t41 = -t66 * t68 + t58;
t35 = t37 ^ 2;
t28 = t71 * t36;
t27 = t68 * t36;
t23 = t68 * t29;
t19 = t20 * t68;
t16 = (-t64 + t65) * t37;
t11 = t12 * t71;
t2 = [1, 0, 0, t70 ^ 2, t70 * t88, 0, 0, 0, pkin(1) * t88, -0.2e1 * pkin(1) * t70, t35, t87, 0, 0, 0, t36 * t119, t37 * t119, t65 * t35, -0.2e1 * t35 * t98, 0.2e1 * t36 * t29, t68 * t87, t36 ^ 2, 0.2e1 * t20 * t99 - 0.2e1 * t36 * t97, -0.2e1 * t10 * t36 + 0.2e1 * t20 * t29 (-t68 * t86 + t7 * t71) * t120, -0.2e1 * t12 * t99 + 0.2e1 * t7 * t36, -0.2e1 * t12 * t29 + 0.2e1 * t36 * t86, t12 ^ 2 + t7 ^ 2 + t86 ^ 2 (t3 * t71 - t5 * t68) * t120, -0.2e1 * t8 * t29 + 0.2e1 * t5 * t36, -0.2e1 * t3 * t36 + 0.2e1 * t8 * t99, t3 ^ 2 + t5 ^ 2 + t8 ^ 2; 0, 0, 0, 0, 0, t70, t105, 0, -t70 * pkin(7), -t85, 0, 0, t37, -t36, 0, -t20, -t21, t23, t16, t27, t28, 0, t68 * t77 - t102, t71 * t77 + t19, t1, t68 * t78 + t11, t71 * t78 - t103, t1 * t55 + t12 * t31 (t33 * t71 - t34 * t68) * t37 + t114, -t22 * t29 + t34 * t36 - t108, t22 * t99 - t33 * t36 - t107, t8 * t22 + t3 * t33 + t5 * t34; 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 1, 0.2e1 * t84, -0.2e1 * t109, t64, t52, 0, 0, 0, t56 * t116, 0.2e1 * t56 * t68, 0.2e1 * t92, t31 * t115, t31 * t117, t90 * t55 ^ 2 + t31 ^ 2, 0.2e1 * t95, t22 * t117, t22 * t116, t22 ^ 2 + t33 ^ 2 + t34 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t37, -t36, 0, -t20, -t21, t23, t16, t27, t28, 0, t68 * t82 - t102, t71 * t82 + t19, t1, t68 * t81 + t11, t71 * t81 - t103, pkin(9) * t1 + t12 * t42 (t44 * t71 - t46 * t68) * t37 + t114, -t30 * t29 + t46 * t36 - t108, t30 * t99 - t44 * t36 - t107, t3 * t44 + t8 * t30 + t5 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, t84, -t109, t64, t52, 0, 0, 0, t106 * t71, -t106 * t68, t91 + t92, t94 * t71, -t94 * t68, pkin(9) * t92 + t31 * t42, t93 + t95, t96 * t68, t96 * t71, t22 * t30 + t33 * t44 + t34 * t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t64, t52, 0, 0, 0, pkin(3) * t115, pkin(3) * t117, 0.2e1 * t91, t42 * t115, t42 * t117, t90 * pkin(9) ^ 2 + t42 ^ 2, 0.2e1 * t93, t30 * t117, t30 * t116, t30 ^ 2 + t44 ^ 2 + t46 ^ 2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, -t99, t36, -t97, -t10, t80 * t37, t97 - 0.2e1 * t111, t86 + t32, -t7 * pkin(4) + qJ(5) * t86, t121 * t37, 0.2e1 * t32 + t76, t36 * t118 + t75, t5 * qJ(5) - t3 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t71, 0, -t49, -t51, t122, t49, t51, t122 * t55, t41, t34, -t33, t34 * qJ(5) - t33 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t71, 0, -t60, -t62, t122, t60, t62, t122 * pkin(9), t41, t46, -t44, t46 * qJ(5) - t44 * t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(4), t72, pkin(4) ^ 2 + t73, 0, t72, t118, t66 ^ 2 + t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t29, t36, 0, t7, t29, 0, -t36, t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, t49, t68, 0, 0, t33; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, 0, 0, t60, t68, 0, 0, t44; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(4), 0, 0, -1, -t66; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, t36, 0, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, 0, 0, t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, qJ(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1;];
MM_reg  = t2;
