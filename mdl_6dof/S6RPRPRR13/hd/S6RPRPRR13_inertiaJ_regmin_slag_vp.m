% Calculate minimal parameter regressor of joint inertia matrix for
% S6RPRPRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d5,d6,theta2]';
% 
% Output:
% MM_reg [((6+1)*6/2)x32]
%   minimal parameter regressor of joint inertia matrix
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MM_reg = S6RPRPRR13_inertiaJ_regmin_slag_vp(qJ, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR13_inertiaJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRPRR13_inertiaJ_regmin_slag_vp: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From inertia_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 20:57:05
% EndTime: 2019-05-05 20:57:09
% DurationCPUTime: 1.17s
% Computational Cost: add. (1582->142), mult. (4382->291), div. (0->0), fcn. (5130->12), ass. (0->96)
t64 = cos(pkin(7));
t61 = sin(pkin(7));
t83 = cos(pkin(6));
t79 = t83 * t61;
t62 = sin(pkin(6));
t63 = cos(pkin(12));
t99 = t62 * t63;
t114 = t64 * t99 + t79;
t60 = sin(pkin(12));
t80 = pkin(1) * t83;
t84 = qJ(2) * t62;
t43 = t60 * t80 + t63 * t84;
t27 = t114 * pkin(9) + t43;
t102 = t60 * t62;
t51 = t63 * t80;
t31 = t83 * pkin(2) + t51 + (-pkin(9) * t64 - qJ(2)) * t102;
t36 = (-pkin(9) * t60 * t61 - pkin(2) * t63 - pkin(1)) * t62;
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t14 = -t67 * t27 + (t31 * t64 + t36 * t61) * t70;
t29 = t67 * t102 - t114 * t70;
t40 = t61 * t99 - t83 * t64;
t66 = sin(qJ(5));
t69 = cos(qJ(5));
t19 = -t29 * t69 - t40 * t66;
t113 = -0.2e1 * t19;
t98 = t64 * t67;
t30 = t67 * t79 + (t60 * t70 + t63 * t98) * t62;
t112 = -0.2e1 * t30;
t68 = cos(qJ(6));
t111 = 0.2e1 * t68;
t110 = 2 * qJ(4);
t109 = pkin(3) + pkin(10);
t55 = t62 ^ 2;
t108 = pkin(1) * t55;
t8 = t30 * pkin(4) + t109 * t40 - t14;
t18 = -t61 * t31 + t64 * t36;
t73 = -t30 * qJ(4) + t18;
t9 = t109 * t29 + t73;
t5 = -t66 * t9 + t69 * t8;
t3 = -t30 * pkin(5) - t5;
t65 = sin(qJ(6));
t107 = t3 * t65;
t106 = t3 * t68;
t105 = t40 * pkin(3);
t20 = t29 * t66 - t40 * t69;
t17 = t20 * t68 + t30 * t65;
t104 = t17 * t65;
t103 = t30 * t66;
t101 = t61 * t67;
t100 = t61 * t70;
t97 = t65 * t19;
t96 = t65 * t66;
t95 = t65 * t68;
t94 = t65 * t69;
t93 = t66 * t109;
t92 = t68 * t19;
t91 = t68 * t66;
t53 = t68 * t69;
t90 = t68 * t109;
t89 = t69 * t17;
t88 = t69 * t19;
t87 = t69 * t66;
t86 = t69 * t109;
t57 = t66 ^ 2;
t59 = t69 ^ 2;
t85 = -t57 - t59;
t39 = qJ(4) * t40;
t82 = -0.2e1 * t87;
t15 = t36 * t101 + t70 * t27 + t31 * t98;
t78 = -t39 + t15;
t77 = -pkin(5) * t69 - pkin(11) * t66;
t6 = t66 * t8 + t69 * t9;
t75 = t40 * t100 - t64 * t29;
t74 = t40 * t101 + t64 * t30;
t10 = -t29 * pkin(4) + t78;
t58 = t68 ^ 2;
t56 = t65 ^ 2;
t47 = t66 * pkin(5) - t69 * pkin(11) + qJ(4);
t45 = -t66 * t100 + t69 * t64;
t44 = t69 * t100 + t66 * t64;
t42 = -t60 * t84 + t51;
t38 = t65 * t47 - t66 * t90;
t37 = t68 * t47 + t65 * t93;
t34 = t65 * t101 + t68 * t45;
t33 = t68 * t101 - t65 * t45;
t28 = t30 ^ 2;
t26 = t30 * t69;
t16 = t20 * t65 - t30 * t68;
t13 = -t14 + t105;
t11 = t29 * pkin(3) + t73;
t7 = t19 * pkin(5) - t20 * pkin(11) + t10;
t4 = t30 * pkin(11) + t6;
t2 = t68 * t4 + t65 * t7;
t1 = -t65 * t4 + t68 * t7;
t12 = [1, 0, 0, 0.2e1 * t63 * t108 + 0.2e1 * t42 * t83, -0.2e1 * t60 * t108 - 0.2e1 * t43 * t83, 0.2e1 * (-t42 * t60 + t43 * t63) * t62, t55 * pkin(1) ^ 2 + t42 ^ 2 + t43 ^ 2, t28, t29 * t112, t40 * t112, 0.2e1 * t29 * t40, t40 ^ 2, -0.2e1 * t14 * t40 + 0.2e1 * t18 * t29, 0.2e1 * t15 * t40 + 0.2e1 * t18 * t30, 0.2e1 * t13 * t30 - 0.2e1 * t29 * t78, -0.2e1 * t11 * t29 - 0.2e1 * t13 * t40, -0.2e1 * t11 * t30 - 0.2e1 * t40 * t78, t11 ^ 2 + t13 ^ 2 + t78 ^ 2, t20 ^ 2, t20 * t113, 0.2e1 * t20 * t30, t19 * t112, t28, 0.2e1 * t10 * t19 + 0.2e1 * t5 * t30, 0.2e1 * t10 * t20 - 0.2e1 * t6 * t30, t17 ^ 2, -0.2e1 * t17 * t16, 0.2e1 * t17 * t19, t16 * t113, t19 ^ 2, 0.2e1 * t1 * t19 + 0.2e1 * t3 * t16, 0.2e1 * t3 * t17 - 0.2e1 * t2 * t19; 0, 0, 0, -t99, t102, 0, -t62 * pkin(1), 0, 0, 0, 0, 0, -t75, t74 (-t29 * t67 - t30 * t70) * t61, t75, -t74, t11 * t64 + (-t13 * t70 + t67 * t78) * t61, 0, 0, 0, 0, 0, t19 * t101 - t44 * t30, t20 * t101 - t45 * t30, 0, 0, 0, 0, 0, t44 * t16 + t33 * t19, t44 * t17 - t34 * t19; 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64 ^ 2 + (t67 ^ 2 + t70 ^ 2) * t61 ^ 2, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t29, -t40, t14, -t15, -t30 * pkin(3) - qJ(4) * t29, -t14 + 0.2e1 * t105, t78 - t39, -t13 * pkin(3) + qJ(4) * t78, t20 * t69, -t20 * t66 - t88, t26, -t103, 0, qJ(4) * t19 + t10 * t66 - t30 * t86, qJ(4) * t20 + t10 * t69 + t30 * t93, t68 * t89 (-t16 * t68 - t104) * t69, t17 * t66 + t68 * t88, -t16 * t66 - t65 * t88, t19 * t66, t1 * t66 + t37 * t19 + (t109 * t16 + t107) * t69, -t38 * t19 - t2 * t66 + (t109 * t17 + t106) * t69; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, -t101, 0, -t100, t101 (pkin(3) * t70 + qJ(4) * t67) * t61, 0, 0, 0, 0, 0, t66 * t101, t69 * t101, 0, 0, 0, 0, 0, t33 * t66 + t44 * t94, -t34 * t66 + t44 * t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, -0.2e1 * pkin(3), t110, pkin(3) ^ 2 + (qJ(4) ^ 2) t59, t82, 0, 0, 0, t66 * t110, t69 * t110, t58 * t59, -0.2e1 * t59 * t95, t87 * t111, t65 * t82, t57, 0.2e1 * t109 * t59 * t65 + 0.2e1 * t37 * t66, -0.2e1 * t38 * t66 + 0.2e1 * t59 * t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t30, -t40, 0, t13, 0, 0, 0, 0, 0, t26, -t103, 0, 0, 0, 0, 0, -t69 * t16 - t19 * t96, -t19 * t91 - t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t100, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, -pkin(3), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t85 * t65, t85 * t68; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t20, -t19, t30, t5, -t6, t104, -t65 * t16 + t17 * t68, t97, t92, 0, -pkin(5) * t16 - pkin(11) * t97 - t106, -pkin(5) * t17 - pkin(11) * t92 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t44, -t45, 0, 0, 0, 0, 0, -t44 * t68, t44 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t66, 0, -t86, t93, t65 * t53 (-t56 + t58) * t69, t96, t91, 0, t65 * t77 - t68 * t86, t65 * t86 + t68 * t77; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t69, -t66, 0, 0, 0, 0, 0, t53, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0, t56, 0.2e1 * t95, 0, 0, 0, pkin(5) * t111, -0.2e1 * pkin(5) * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t17, -t16, t19, t1, -t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, -t34; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, -t94, t66, t37, -t38; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t96, -t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t68, 0, -t65 * pkin(11), -t68 * pkin(11); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 1, 0, 0;];
MM_reg  = t12;
