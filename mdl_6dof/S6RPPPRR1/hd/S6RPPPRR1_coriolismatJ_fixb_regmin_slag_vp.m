% Calculate minimal parameter regressor of coriolis matrix for
% S6RPPPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2]';
% 
% Output:
% cmat_reg [(6*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S6RPPPRR1_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR1_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:30:17
% EndTime: 2019-03-09 01:30:19
% DurationCPUTime: 0.94s
% Computational Cost: add. (414->142), mult. (918->191), div. (0->0), fcn. (714->6), ass. (0->115)
t69 = sin(qJ(6));
t136 = 0.2e1 * t69;
t63 = t69 ^ 2;
t71 = cos(qJ(6));
t65 = t71 ^ 2;
t43 = t65 - t63;
t72 = cos(qJ(5));
t126 = t71 * t72;
t90 = t126 * t136;
t74 = qJD(1) * t90 - t43 * qJD(5);
t70 = sin(qJ(5));
t64 = t70 ^ 2;
t135 = t64 / 0.2e1;
t134 = t70 * pkin(8);
t133 = t72 * pkin(5);
t56 = sin(pkin(9)) * pkin(1) + qJ(3);
t40 = -pkin(7) + t56;
t132 = t40 * t69;
t66 = t72 ^ 2;
t131 = t66 * t69;
t130 = t66 * t71;
t32 = t133 + t134;
t129 = t69 * t32;
t128 = t71 * t32;
t127 = t71 * t40;
t104 = t72 * qJD(5);
t52 = t69 * t104;
t114 = qJD(6) * t71;
t55 = t70 * t114;
t125 = t55 + t52;
t42 = t64 - t66;
t124 = t64 + t66;
t101 = t72 * t132;
t41 = cos(pkin(9)) * pkin(1) + pkin(2) + qJ(4);
t85 = t70 * pkin(5) - t72 * pkin(8);
t73 = t41 + t85;
t7 = t70 * t132 - t71 * t73;
t1 = -t7 * t72 + (t101 + t128) * t70;
t123 = t1 * qJD(1);
t100 = t40 * t126;
t102 = t70 * t127;
t8 = t69 * t73 + t102;
t2 = t8 * t72 + (-t100 + t129) * t70;
t122 = t2 * qJD(1);
t3 = -t40 * t131 - t7 * t70;
t121 = t3 * qJD(1);
t4 = -t66 * t127 - t8 * t70;
t120 = t4 * qJD(1);
t119 = qJD(3) * t70;
t118 = qJD(3) * t72;
t117 = qJD(4) * t70;
t116 = qJD(5) * t71;
t115 = qJD(6) * t69;
t88 = 0.1e1 / 0.2e1 + t135 + t66 / 0.2e1;
t12 = t88 * t69;
t113 = t12 * qJD(1);
t13 = t88 * t71;
t112 = t13 * qJD(1);
t28 = t124 * t69;
t111 = t28 * qJD(1);
t29 = t42 * t69;
t110 = t29 * qJD(1);
t30 = t124 * t71;
t109 = t30 * qJD(1);
t31 = t42 * t71;
t108 = t31 * qJD(1);
t107 = t41 * qJD(1);
t106 = t42 * qJD(1);
t105 = t56 * qJD(1);
t59 = t70 * qJD(1);
t58 = t70 * qJD(5);
t60 = t72 * qJD(1);
t103 = t72 * qJD(6);
t99 = t69 * t103;
t98 = t71 * t103;
t97 = t69 * t114;
t96 = t69 * t116;
t95 = t69 * t60;
t94 = t71 * t104;
t93 = t71 * t60;
t92 = t70 * t104;
t91 = t70 * t60;
t89 = -qJD(3) + t107;
t86 = qJD(5) * t90;
t54 = t70 * t115;
t84 = t54 - t94;
t83 = (-qJD(6) - t59) * t72;
t82 = t134 / 0.2e1 + t133 / 0.2e1;
t77 = t32 / 0.2e1 + t82;
t9 = t77 * t69;
t81 = pkin(5) * t116 - t9 * qJD(1);
t10 = t77 * t71;
t80 = pkin(5) * t69 * qJD(5) + t10 * qJD(1);
t79 = t71 * t83;
t20 = (t63 / 0.2e1 - t65 / 0.2e1) * t72;
t78 = -t20 * qJD(1) + t96;
t76 = t69 * qJD(1) * t130 + t20 * qJD(5);
t27 = t43 * t66;
t75 = t27 * qJD(1) + t86;
t57 = t104 / 0.2e1;
t53 = t71 * t59;
t51 = t69 * t58;
t50 = t69 * t59;
t33 = t56 * qJD(3);
t26 = t53 + t114;
t25 = t50 + t115;
t23 = (t59 + qJD(6) / 0.2e1) * t72;
t19 = t51 - t98;
t18 = t71 * t58 + t99;
t17 = t20 * qJD(6);
t15 = -t130 / 0.2e1 + (-t64 / 0.2e1 + 0.1e1 / 0.2e1) * t71;
t14 = t131 / 0.2e1 + (t135 - 0.1e1 / 0.2e1) * t69;
t6 = -t101 + t128 / 0.2e1 - t82 * t71;
t5 = -t100 - t129 / 0.2e1 + t82 * t69;
t11 = [0, 0, 0, 0, 0, qJD(3), t33, qJD(3), qJD(4), t41 * qJD(4) + t33, -t92, t42 * qJD(5), 0, 0, 0, t41 * t104 + t117, qJD(4) * t72 - t41 * t58, -t65 * t92 - t66 * t97, -t27 * qJD(6) + t70 * t86, -t31 * qJD(5) - t70 * t99, t29 * qJD(5) - t70 * t98, t92, -t28 * qJD(3) + t1 * qJD(5) + t4 * qJD(6) + t71 * t117, -t30 * qJD(3) - t2 * qJD(5) - t3 * qJD(6) - t69 * t117; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, qJD(1), t105, qJD(1), 0, t105, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t109; 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), t107, 0, 0, 0, 0, 0, t59, t60, 0, 0, 0, 0, 0, t15 * qJD(6) + t53, t14 * qJD(6) - t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t91, t106, -t58, -t104, 0, -t40 * t58 + t41 * t60, -t40 * t104 - t41 * t59, -t17 + (-t65 * t60 - t96) * t70, t70 * t74 - 0.2e1 * t72 * t97, t52 - t108, t94 + t110, t23, t123 + (t85 * t69 - t102) * qJD(5) + t6 * qJD(6), -t122 + (-pkin(8) * t126 + (pkin(5) * t71 + t132) * t70) * qJD(5) + t5 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, -t75, t69 * t83, t79, t57, t15 * qJD(4) + t6 * qJD(5) - t8 * qJD(6) + t120, t14 * qJD(4) + t5 * qJD(5) + t7 * qJD(6) - t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, t58, 0, 0, 0, 0, 0, t84, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t19, t18; 0, 0, 0, 0, 0, -qJD(1), -t105, -qJD(1), 0, -t105 - qJD(4), 0, 0, 0, 0, 0, -t104, t58, 0, 0, 0, 0, 0, t84 + t111, t109 + t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t60, t59, 0, 0, 0, 0, 0, -t93, t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t25, t26; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t89, 0, 0, 0, 0, 0, -t59, -t60, 0, 0, 0, 0, 0, -t13 * qJD(6) - t53, t12 * qJD(6) + t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t104, 0, 0, 0, 0, 0, -t18, t19; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 - t125, t84 + t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t91, -t106, 0, 0, 0, -t89 * t72, t89 * t70, t65 * t91 - t17, t79 * t136, t55 + t108, -t54 - t110, -t23, -t10 * qJD(6) + t71 * t118 - t123, t9 * qJD(6) - t69 * t118 + t122; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t60, -t59, 0, 0, 0, 0, 0, t93, -t95; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t97, t43 * qJD(6), 0, 0, 0, -pkin(5) * t115, -pkin(5) * t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, -t74, t26, -t25, -t60 / 0.2e1, -pkin(8) * t114 - t80, pkin(8) * t115 - t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, t75 (t95 - t116) * t70, t71 * t91 + t51, t57, t13 * qJD(4) + t10 * qJD(5) - t69 * t119 - t120, -t12 * qJD(4) - t9 * qJD(5) - t71 * t119 + t121; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t50, -t53; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t112, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t78, t74, -t53, t50, t60 / 0.2e1, t80, t81; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg  = t11;
