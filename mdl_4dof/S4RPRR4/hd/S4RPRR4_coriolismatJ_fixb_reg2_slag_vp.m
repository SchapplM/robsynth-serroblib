% Calculate inertial parameters regressor of coriolis matrix for
% S4RPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
% 
% Output:
% cmat_reg [(4*4)x(4*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:50:38
% EndTime: 2019-12-31 16:50:40
% DurationCPUTime: 0.83s
% Computational Cost: add. (873->122), mult. (1851->197), div. (0->0), fcn. (1518->6), ass. (0->106)
t75 = sin(qJ(4));
t141 = 0.2e1 * t75;
t76 = sin(qJ(3));
t140 = t76 * pkin(3);
t78 = cos(qJ(3));
t139 = t78 * pkin(6);
t66 = sin(pkin(7)) * pkin(1) + pkin(5);
t134 = t66 * t78;
t106 = t75 * t134;
t77 = cos(qJ(4));
t67 = -cos(pkin(7)) * pkin(1) - pkin(2);
t91 = -t78 * pkin(3) - t76 * pkin(6);
t79 = t67 + t91;
t29 = -t77 * t79 + t106;
t138 = t29 * t78;
t128 = t77 * t66;
t104 = t78 * t128;
t30 = t75 * t79 + t104;
t137 = t30 * t78;
t55 = -t139 + t140;
t129 = t77 * t55;
t131 = t75 * t76;
t52 = t66 * t131;
t35 = t52 + t129;
t136 = t35 * t76;
t105 = t76 * t128;
t132 = t75 * t55;
t36 = -t105 + t132;
t135 = t36 * t76;
t71 = t75 ^ 2;
t133 = t71 * t78;
t130 = t76 * t78;
t72 = t76 ^ 2;
t127 = t77 * t72;
t73 = t77 ^ 2;
t126 = t71 + t73;
t60 = t73 - t71;
t74 = t78 ^ 2;
t61 = t74 - t72;
t3 = (t136 - t138) * t77 + (t135 + t137) * t75;
t125 = t3 * qJD(1);
t6 = t29 * t76 + (t35 - 0.2e1 * t52) * t78;
t124 = t6 * qJD(1);
t7 = t36 * t78 + (-t30 + 0.2e1 * t104) * t76;
t123 = t7 * qJD(1);
t122 = qJD(1) * t78;
t121 = qJD(3) * t75;
t120 = qJD(3) * t77;
t119 = qJD(4) * t75;
t118 = qJD(4) * t77;
t117 = qJD(4) * t78;
t14 = -t72 * t66 * t75 - t138;
t116 = t14 * qJD(1);
t15 = -t66 * t127 - t137;
t115 = t15 * qJD(1);
t44 = (t71 / 0.2e1 - t73 / 0.2e1) * t76;
t114 = t44 * qJD(4);
t49 = t61 * t75;
t113 = t49 * qJD(1);
t50 = t77 * t74 - t127;
t112 = t50 * qJD(1);
t111 = t61 * qJD(1);
t110 = t76 * qJD(1);
t109 = t76 * qJD(3);
t108 = t76 * qJD(4);
t107 = t78 * qJD(3);
t103 = t77 * t110;
t102 = t75 * t117;
t101 = t77 * t117;
t100 = t67 * t110;
t99 = t67 * t122;
t98 = t75 * t118;
t97 = t75 * t120;
t96 = t76 * t107;
t65 = t78 * t110;
t95 = t77 * t109;
t94 = t72 * t98;
t93 = t75 * t95;
t92 = -qJD(4) + t122;
t90 = -t35 * t75 + t36 * t77;
t2 = (t135 / 0.2e1 + t137 / 0.2e1) * t77 + (-t136 / 0.2e1 + t138 / 0.2e1) * t75 + (t72 / 0.2e1 - t74 / 0.2e1) * t66;
t4 = t66 ^ 2 * t130 - t29 * t35 + t30 * t36;
t89 = t4 * qJD(1) + t2 * qJD(2);
t39 = (-0.1e1 + t126) * t130;
t88 = -t2 * qJD(1) - t39 * qJD(2);
t87 = t92 * t76;
t86 = t139 / 0.2e1 - t140 / 0.2e1;
t82 = -t55 / 0.2e1 + t86;
t31 = t82 * t75;
t85 = pkin(3) * t120 + t31 * qJD(1);
t32 = t82 * t77;
t84 = pkin(3) * t121 - t32 * qJD(1);
t83 = t77 * t87;
t38 = -t44 * qJD(1) + t97;
t33 = t75 * qJD(1) * t127 + t44 * qJD(3);
t48 = t60 * t72;
t81 = t48 * qJD(1) + 0.2e1 * t93;
t80 = -t60 * qJD(3) + t103 * t141;
t69 = t73 * t78;
t68 = t109 / 0.2e1;
t64 = t75 * t109;
t47 = t65 - t108 / 0.2e1;
t17 = t52 + t129 / 0.2e1 + t86 * t77;
t16 = t105 - t132 / 0.2e1 - t86 * t75;
t1 = t2 * qJD(3);
t5 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t96, t61 * qJD(3), 0, -t96, 0, 0, t67 * t109, t67 * t107, 0, 0, t73 * t96 - t94, -t48 * qJD(4) - 0.2e1 * t78 * t93, -t50 * qJD(3) + t76 * t102, t71 * t96 + t94, t49 * qJD(3) + t76 * t101, -t96, -t6 * qJD(3) - t15 * qJD(4), t7 * qJD(3) + t14 * qJD(4), -t3 * qJD(3), t4 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t65, t111, t107, -t65, -t109, 0, -t66 * t107 + t100, t66 * t109 + t99, 0, 0, -t114 + (t73 * t110 + t97) * t78, (t69 - t133) * qJD(3) + 0.2e1 * (-qJD(4) - t122) * t77 * t131, t64 - t112, t114 + (t71 * t110 - t97) * t78, t95 + t113, -t47, -t124 + (t91 * t75 - t104) * qJD(3) + t17 * qJD(4), t123 + (t91 * t77 + t106) * qJD(3) + t16 * qJD(4), t90 * qJD(3) - t125, (-pkin(3) * t134 + t90 * pkin(6)) * qJD(3) + t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t33, -t81, t75 * t87, t33, t83, t68, t17 * qJD(3) - t30 * qJD(4) - t115, t16 * qJD(3) + t29 * qJD(4) + t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t39 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t109, -t107, 0, 0, 0, 0, 0, 0, 0, 0, -t95 - t102, t64 - t101, (t69 + t133) * qJD(3), (t126 * t139 - t140) * qJD(3) - t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t75 * t107 - t77 * t108, -t77 * t107 + t75 * t108, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t65, -t111, 0, t65, 0, 0, -t100, -t99, 0, 0, -t73 * t65 - t114, t83 * t141, -t101 + t112, -t71 * t65 + t114, t102 - t113, t47, t32 * qJD(4) + t124, -t31 * qJD(4) - t123, t125, -t89; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t98, t60 * qJD(4), 0, -t98, 0, 0, -pkin(3) * t119, -pkin(3) * t118, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, -t80, -t92 * t77, -t38, t92 * t75, -t110 / 0.2e1, -pkin(6) * t118 - t84, pkin(6) * t119 - t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t81, (-t75 * t110 + t120) * t78, -t33, (-t103 - t121) * t78, t68, -t32 * qJD(3) + t115, t31 * qJD(3) - t116, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t38, t80, t77 * t122, t38, -t75 * t122, t110 / 0.2e1, t84, t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t5;
