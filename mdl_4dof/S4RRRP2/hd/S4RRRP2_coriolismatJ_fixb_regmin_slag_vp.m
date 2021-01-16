% Calculate minimal parameter regressor of coriolis matrix for
% S4RRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3]';
% 
% Output:
% cmat_reg [(4*%NQJ)%x17]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-15 11:04
% Revision: d12c3222fdeb2c5f3b3c8fa5751e113be2fc3aae (2021-01-15)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S4RRRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-15 11:04:31
% EndTime: 2021-01-15 11:04:32
% DurationCPUTime: 0.54s
% Computational Cost: add. (663->123), mult. (1373->161), div. (0->0), fcn. (951->4), ass. (0->110)
t89 = cos(qJ(2));
t139 = t89 * pkin(1);
t75 = -pkin(2) - t139;
t145 = t75 / 0.2e1 - pkin(2) / 0.2e1;
t114 = qJD(1) + qJD(2);
t86 = sin(qJ(3));
t84 = t86 ^ 2;
t88 = cos(qJ(3));
t85 = t88 ^ 2;
t71 = t85 + t84;
t147 = t114 * t71;
t72 = t85 - t84;
t146 = t114 * t72;
t142 = pkin(3) * t86;
t87 = sin(qJ(2));
t141 = t87 * pkin(1);
t140 = t88 * pkin(3);
t74 = pkin(6) + t141;
t124 = qJ(4) + t74;
t43 = t124 * t86;
t138 = t43 * t86;
t44 = t124 * t88;
t137 = t44 * t88;
t76 = -pkin(2) - t140;
t59 = t76 - t139;
t50 = t59 * t86;
t136 = t59 * t88;
t132 = pkin(6) + qJ(4);
t60 = t132 * t86;
t135 = t60 * t86;
t61 = t132 * t88;
t134 = t61 * t88;
t65 = t76 * t86;
t133 = t76 * t88;
t42 = t71 * t139;
t62 = t71 * qJD(4);
t131 = t42 * qJD(2) + t62;
t129 = pkin(1) * qJD(1);
t109 = t87 * t129;
t68 = t86 * t109;
t128 = pkin(1) * qJD(2);
t112 = t87 * t128;
t70 = t86 * t112;
t130 = t68 + t70;
t127 = pkin(2) * qJD(2);
t3 = pkin(3) * t50;
t126 = t3 * qJD(1);
t15 = t137 + t138;
t8 = (t15 * t89 + t59 * t87) * pkin(1);
t125 = t8 * qJD(1);
t123 = qJD(1) * t75;
t122 = t15 * qJD(1);
t113 = t86 * t140;
t26 = -t50 + t113;
t121 = t26 * qJD(1);
t83 = t84 * pkin(3);
t35 = t83 + t136;
t120 = t35 * qJD(1);
t119 = t42 * qJD(1);
t118 = t44 * qJD(3);
t117 = t61 * qJD(3);
t81 = t86 * qJD(3);
t116 = t86 * qJD(4);
t82 = t88 * qJD(3);
t115 = t88 * qJD(4);
t111 = pkin(3) * t82;
t110 = pkin(3) * t116;
t107 = t141 / 0.2e1;
t106 = -t139 / 0.2e1;
t105 = t86 * t123;
t104 = t88 * t123;
t103 = -t65 / 0.2e1 - t50 / 0.2e1;
t102 = pkin(1) * t114;
t101 = t88 * t112;
t100 = t86 * t106;
t99 = t88 * t106;
t57 = t114 * t86;
t58 = t114 * t88;
t17 = t134 + t135;
t92 = t106 - t76 / 0.2e1 - t59 / 0.2e1;
t1 = t92 * t142;
t9 = pkin(3) * t65;
t98 = -t1 * qJD(1) + t9 * qJD(2);
t4 = t107 + (-t61 / 0.2e1 - t44 / 0.2e1) * t88 + (-t60 / 0.2e1 - t43 / 0.2e1) * t86;
t97 = -t4 * qJD(1) + t17 * qJD(2);
t69 = t88 * t109;
t96 = -t69 - t101;
t10 = (t106 + t140) * t86 + t103;
t36 = -t65 + t113;
t95 = -t10 * qJD(1) - t36 * qJD(2);
t12 = t92 * t88 - t83;
t49 = t83 + t133;
t94 = -t12 * qJD(1) + t49 * qJD(2);
t93 = t106 - t145;
t22 = t93 * t86;
t91 = t22 * qJD(1) + t86 * t127;
t23 = t93 * t88;
t90 = t23 * qJD(1) + t88 * t127;
t80 = pkin(3) * t81;
t73 = t86 * t82;
t63 = t72 * qJD(3);
t51 = pkin(3) * t57;
t39 = t86 * t58;
t25 = t145 * t88 + t99;
t24 = t145 * t86 + t100;
t13 = t83 + t133 / 0.2e1 + t136 / 0.2e1 + t99;
t11 = (t106 - t140) * t86 - t103;
t5 = t134 / 0.2e1 + t137 / 0.2e1 + t135 / 0.2e1 + t138 / 0.2e1 + t107;
t2 = pkin(3) * t100 + (t59 + t76) * t142 / 0.2e1;
t6 = [0, 0, 0, 0, -t112, -t89 * t128, t73, t63, 0, 0, 0, t75 * t81 - t101, t75 * t82 + t70, -t26 * qJD(3) - t101, t35 * qJD(3) + t70, t131, t8 * qJD(2) + t3 * qJD(3) + t15 * qJD(4); 0, 0, 0, 0, -t87 * t102, -t89 * t102, t73, t63, 0, 0, 0, t24 * qJD(3) + t96, t25 * qJD(3) + t130, t11 * qJD(3) + t96, t13 * qJD(3) + t130, t119 + t131, t125 + t2 * qJD(3) + t5 * qJD(4) + (t17 * t89 + t76 * t87) * t128; 0, 0, 0, 0, 0, 0, t39, t146, t82, -t81, 0, t24 * qJD(2) - t74 * t82 + t105, t25 * qJD(2) + t74 * t81 + t104, t11 * qJD(2) - t118 - t121, t13 * qJD(2) + t43 * qJD(3) + t120, -t111, -pkin(3) * t118 + t2 * qJD(2) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t5 * qJD(2) + t122; 0, 0, 0, 0, t109, t89 * t129, t73, t63, 0, 0, 0, -t22 * qJD(3) + t69, -t23 * qJD(3) - t68, -t10 * qJD(3) + t69, -t12 * qJD(3) - t68, t62 - t119, -t1 * qJD(3) - t4 * qJD(4) - t125; 0, 0, 0, 0, 0, 0, t73, t63, 0, 0, 0, -pkin(2) * t81, -pkin(2) * t82, -t36 * qJD(3), t49 * qJD(3), t62, t9 * qJD(3) + t17 * qJD(4); 0, 0, 0, 0, 0, 0, t39, t146, t82, -t81, 0, -pkin(6) * t82 - t91, pkin(6) * t81 - t90, t95 - t117, t60 * qJD(3) + t94, -t111, -pkin(3) * t117 + t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t147, t97; 0, 0, 0, 0, 0, 0, -t39, -t146, 0, 0, 0, t22 * qJD(2) - t105, t23 * qJD(2) - t104, t10 * qJD(2) - t116 + t121, t12 * qJD(2) - t115 - t120, 0, t1 * qJD(2) - t110 - t126; 0, 0, 0, 0, 0, 0, -t39, -t146, 0, 0, 0, t91, t90, -t95 - t116, -t94 - t115, 0, -t98 - t110; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t57, -t58, 0, -t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t82, -t147, t4 * qJD(2) - t122 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t81, t82, -t147, t80 - t97; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t57, t58, 0, t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t6;
