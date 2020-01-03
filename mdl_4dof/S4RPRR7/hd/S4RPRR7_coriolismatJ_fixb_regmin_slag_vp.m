% Calculate minimal parameter regressor of coriolis matrix for
% S4RPRR7
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
% cmat_reg [(4*%NQJ)%x21]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RPRR7_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:54:06
% EndTime: 2019-12-31 16:54:08
% DurationCPUTime: 0.67s
% Computational Cost: add. (815->108), mult. (1870->175), div. (0->0), fcn. (1987->6), ass. (0->105)
t133 = cos(qJ(3));
t73 = sin(pkin(7));
t74 = cos(pkin(7));
t76 = sin(qJ(3));
t58 = -t133 * t74 + t76 * t73;
t55 = t58 ^ 2;
t129 = t76 * t74;
t97 = t133 * t73;
t60 = t97 + t129;
t56 = t60 ^ 2;
t137 = t55 + t56;
t104 = t56 - t55;
t75 = sin(qJ(4));
t71 = t75 ^ 2;
t77 = cos(qJ(4));
t72 = t77 ^ 2;
t66 = t72 - t71;
t130 = t75 * t77;
t95 = 0.2e1 * t60 * t130;
t79 = qJD(1) * t95 - t66 * qJD(3);
t134 = t60 * pkin(3);
t135 = t58 * pkin(6);
t36 = t134 + t135;
t136 = t36 / 0.2e1;
t127 = pkin(5) + qJ(2);
t62 = t127 * t74;
t96 = t127 * t73;
t37 = t133 * t96 + t76 * t62;
t132 = t37 * t77;
t38 = t133 * t62 - t76 * t96;
t131 = t75 * t38;
t22 = t75 * t58;
t24 = t75 * t60;
t128 = t77 * t38;
t27 = t77 * t58;
t63 = t73 ^ 2 + t74 ^ 2;
t68 = -t74 * pkin(2) - pkin(1);
t90 = t58 * pkin(3) - t60 * pkin(6);
t78 = t68 + t90;
t9 = -t77 * t78 + t131;
t1 = t36 * t27 + (-t9 + t131) * t60;
t126 = t1 * qJD(1);
t10 = t75 * t78 + t128;
t2 = -t36 * t22 + (-t10 + t128) * t60;
t125 = t2 * qJD(1);
t7 = -t37 * t24 + t9 * t58;
t124 = t7 * qJD(1);
t8 = -t10 * t58 + t60 * t132;
t123 = t8 * qJD(1);
t122 = qJD(2) * t77;
t121 = qJD(3) * t77;
t120 = qJD(4) * t75;
t119 = qJD(4) * t77;
t11 = t104 * t75;
t118 = t11 * qJD(1);
t12 = t137 * t75;
t117 = t12 * qJD(1);
t13 = t104 * t77;
t116 = t13 * qJD(1);
t115 = t104 * qJD(1);
t114 = t22 * qJD(1);
t113 = t24 * qJD(1);
t112 = t27 * qJD(1);
t31 = t137 * t77;
t111 = t31 * qJD(1);
t54 = t97 / 0.2e1 + t129 / 0.2e1;
t110 = t54 * qJD(1);
t109 = t58 * qJD(1);
t53 = t58 * qJD(3);
t108 = t60 * qJD(1);
t107 = t60 * qJD(3);
t61 = t63 * qJ(2);
t106 = t61 * qJD(1);
t105 = t63 * qJD(1);
t103 = t58 * t119;
t102 = t58 * t108;
t101 = t58 * t107;
t100 = t75 * t119;
t99 = t75 * t121;
t98 = t77 * t108;
t94 = qJD(1) * t68 + qJD(2);
t93 = -qJD(4) - t109;
t91 = qJD(3) * t95;
t89 = t93 * t77;
t88 = t135 / 0.2e1 + t134 / 0.2e1;
t82 = t136 + t88;
t3 = t82 * t75;
t87 = pkin(3) * t121 - t3 * qJD(1);
t5 = t82 * t77;
t86 = pkin(3) * t75 * qJD(3) + t5 * qJD(1);
t85 = t60 * t89;
t21 = (t71 / 0.2e1 - t72 / 0.2e1) * t60;
t84 = -t21 * qJD(1) + t99;
t83 = t54 * qJD(4) + t102;
t81 = t56 * qJD(1) * t130 + t21 * qJD(3);
t30 = t66 * t56;
t80 = t30 * qJD(1) + t91;
t52 = t54 * qJD(3);
t51 = t77 * t107;
t18 = t22 * qJD(4);
t17 = t21 * qJD(4);
t14 = -t114 - t120;
t6 = t37 * t75 + (t136 - t88) * t77;
t4 = t132 + (-t36 / 0.2e1 + t88) * t75;
t15 = [0, 0, 0, 0, 0, t63 * qJD(2), t61 * qJD(2), -t101, -t104 * qJD(3), 0, 0, 0, t68 * t107, -t68 * t53, -t56 * t100 - t72 * t101, -t30 * qJD(4) + t58 * t91, -t60 * t58 * t120 + t13 * qJD(3), -t11 * qJD(3) - t60 * t103, t101, t12 * qJD(2) + t1 * qJD(3) + t8 * qJD(4), t31 * qJD(2) + t2 * qJD(3) + t7 * qJD(4); 0, 0, 0, 0, 0, t105, t106, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t117, t111; 0, 0, 0, 0, 0, 0, 0, -t102, -t115, -t53, -t107, 0, -t38 * qJD(3) + t68 * t108, t37 * qJD(3) - t68 * t109, -t17 + (-t72 * t108 - t99) * t58, -0.2e1 * t60 * t100 + t79 * t58, t75 * t107 + t116, t51 - t118, t83, t126 + (t90 * t75 - t128) * qJD(3) + t6 * qJD(4), t125 + (t90 * t77 + t131) * qJD(3) + t4 * qJD(4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t81, -t80, t93 * t24, t85, t52, t6 * qJD(3) - t10 * qJD(4) + t123, t4 * qJD(3) + t9 * qJD(4) + t124; 0, 0, 0, 0, 0, -t105, -t106, 0, 0, 0, 0, 0, t107, -t53, 0, 0, 0, 0, 0, -t18 + t51 - t117, -t24 * qJD(3) - t103 - t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t108, -t109, 0, 0, 0, 0, 0, t98, -t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t14, t89; 0, 0, 0, 0, 0, 0, 0, t102, t115, 0, 0, 0, -t94 * t60, t94 * t58, t72 * t102 - t17, 0.2e1 * t75 * t85, t27 * qJD(4) - t116, -t18 + t118, -t83, -t5 * qJD(4) - t60 * t122 - t126, t24 * qJD(2) + t3 * qJD(4) - t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, t109, 0, 0, 0, 0, 0, -t98, t113; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t100, t66 * qJD(4), 0, 0, 0, -pkin(3) * t120, -pkin(3) * t119; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t79, t112 + t119, t14, -t110, -pkin(6) * t119 - t86, pkin(6) * t120 - t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t80, -t27 * qJD(3) + t75 * t102, t22 * qJD(3) + t58 * t98, t52, t22 * qJD(2) + t5 * qJD(3) - t123, -t3 * qJD(3) + t58 * t122 - t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t114, t77 * t109; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t79, -t112, t114, t110, t86, t87; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t15;
