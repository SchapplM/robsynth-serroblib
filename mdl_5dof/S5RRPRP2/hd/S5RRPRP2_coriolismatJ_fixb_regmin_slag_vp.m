% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x18]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP2_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:45
% EndTime: 2019-12-31 19:49:47
% DurationCPUTime: 0.78s
% Computational Cost: add. (807->130), mult. (1607->156), div. (0->0), fcn. (1271->6), ass. (0->117)
t70 = sin(pkin(8));
t73 = sin(qJ(2));
t61 = t70 * t73 * pkin(1);
t75 = cos(qJ(2));
t129 = t75 * pkin(1);
t64 = pkin(2) + t129;
t71 = cos(pkin(8));
t89 = t71 * t64 - t61;
t44 = -pkin(3) - t89;
t131 = t71 * pkin(2);
t63 = -pkin(3) - t131;
t136 = t63 / 0.2e1 + t44 / 0.2e1;
t72 = sin(qJ(4));
t139 = t72 * t136;
t74 = cos(qJ(4));
t138 = t74 * t136;
t85 = -t74 * pkin(4) - t72 * qJ(5);
t81 = -pkin(3) + t85;
t26 = t81 - t89;
t25 = t26 * t72;
t51 = t81 - t131;
t46 = t51 * t72;
t119 = t25 / 0.2e1 + t46 / 0.2e1;
t113 = t74 * qJ(5);
t130 = t72 * pkin(4);
t55 = t113 - t130;
t123 = t55 * t74;
t137 = -t119 - t123;
t101 = qJD(1) + qJD(2);
t68 = t72 ^ 2;
t69 = t74 ^ 2;
t58 = t69 - t68;
t134 = t101 * t58;
t128 = t26 * t55;
t127 = t26 * t74;
t126 = t51 * t74;
t125 = t55 * t51;
t124 = t55 * t72;
t122 = t71 * t73;
t50 = t71 * t129 - t61;
t121 = t72 * t50;
t120 = t74 * t50;
t110 = qJD(2) * t72;
t49 = (t70 * t75 + t122) * pkin(1);
t36 = t49 * t110;
t65 = t68 * qJD(5);
t117 = t65 - t36;
t116 = pkin(1) * qJD(1);
t115 = pkin(1) * qJD(2);
t19 = (t68 + t69) * t50;
t87 = pkin(1) * t122 + t70 * t64;
t45 = pkin(7) + t87;
t1 = t45 * t19 + t26 * t49;
t114 = t1 * qJD(1);
t112 = qJD(1) * t72;
t111 = qJD(1) * t74;
t109 = qJD(2) * t74;
t10 = -t89 * t49 + t87 * t50;
t108 = t10 * qJD(1);
t15 = -t124 + t127;
t107 = t15 * qJD(1);
t16 = -t25 - t123;
t106 = t16 * qJD(1);
t105 = t19 * qJD(1);
t104 = t72 * qJD(4);
t66 = t74 * qJD(4);
t103 = t74 * qJD(5);
t102 = qJD(4) * qJ(5);
t100 = qJD(1) * t128;
t99 = t26 * t112;
t98 = t44 * t112;
t97 = t44 * t111;
t96 = t49 * t109;
t95 = t45 * t104;
t62 = t70 * pkin(2) + pkin(7);
t94 = t62 * t104;
t93 = t45 * t66;
t92 = t62 * t66;
t91 = t51 / 0.2e1 + t26 / 0.2e1;
t88 = pkin(1) * t101;
t59 = t72 * t103;
t86 = t59 - t96;
t20 = -t124 + t126;
t34 = t120 / 0.2e1;
t5 = t91 * t74 - t124 + t34;
t84 = t5 * qJD(1) + t20 * qJD(2);
t21 = -t46 - t123;
t33 = -t121 / 0.2e1;
t4 = t33 + t137;
t83 = t4 * qJD(1) + t21 * qJD(2);
t82 = t55 * qJD(4) + t72 * qJD(5);
t76 = (t113 / 0.2e1 - t130 / 0.2e1) * t50;
t2 = t91 * t55 + t76;
t80 = t2 * qJD(1) + qJD(2) * t125;
t32 = t121 / 0.2e1;
t8 = t32 + t119;
t79 = t8 * qJD(1) + t51 * t110;
t11 = t33 - t139;
t78 = t11 * qJD(1) - t63 * t110;
t35 = -t120 / 0.2e1;
t12 = t35 - t138;
t77 = t12 * qJD(1) - t63 * t109;
t47 = t85 * qJD(4) + t103;
t60 = t72 * t66;
t57 = t58 * qJD(4);
t54 = t101 * t68;
t48 = t101 * t74 * t72;
t38 = t49 * t111;
t37 = t49 * t112;
t17 = t19 * qJD(2);
t14 = t35 + t138;
t13 = t33 + t139;
t9 = t32 - t119;
t7 = t33 - t137;
t6 = -t127 / 0.2e1 + t124 + t34 - t126 / 0.2e1;
t3 = -t125 / 0.2e1 - t128 / 0.2e1 + t76;
t18 = [0, 0, 0, 0, -t73 * t115, -t75 * t115, t10 * qJD(2), t60, t57, 0, 0, 0, t44 * t104 - t96, t44 * t66 + t36, -t16 * qJD(4) + t86, t17, -t15 * qJD(4) + t117, t1 * qJD(2) - t82 * t26; 0, 0, 0, 0, -t73 * t88, -t75 * t88, t108 + (-t49 * t71 + t50 * t70) * qJD(2) * pkin(2), t60, t57, 0, 0, 0, t13 * qJD(4) - t38 - t96, t14 * qJD(4) + t36 + t37, t7 * qJD(4) - t38 + t86, t17 + t105, t6 * qJD(4) + t117 - t37, t114 + (t62 * t19 + t49 * t51) * qJD(2) + t3 * qJD(4) + t9 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t48, t134, t66, -t104, 0, t13 * qJD(2) - t93 + t98, t14 * qJD(2) + t95 + t97, t7 * qJD(2) - t106 - t93, t47, t6 * qJD(2) - t107 - t95, t3 * qJD(2) + t47 * t45 - t100; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t66, t54, t9 * qJD(2) + t93 - t99; 0, 0, 0, 0, t73 * t116, t75 * t116, -t108, t60, t57, 0, 0, 0, -t11 * qJD(4) + t38, -t12 * qJD(4) - t37, -t4 * qJD(4) + t38 + t59, -t105, -t5 * qJD(4) + t37 + t65, -t2 * qJD(4) - t8 * qJD(5) - t114; 0, 0, 0, 0, 0, 0, 0, t60, t57, 0, 0, 0, t63 * t104, t63 * t66, -t21 * qJD(4) + t59, 0, -t20 * qJD(4) + t65, -t82 * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, t48, t134, t66, -t104, 0, -t78 - t92, -t77 + t94, -t83 - t92, t47, -t84 - t94, t47 * t62 - t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, t66, t54, -t79 + t92; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104, -t66, -t104, 0, t66, t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104; 0, 0, 0, 0, 0, 0, 0, -t48, -t134, 0, 0, 0, t11 * qJD(2) - t98, t12 * qJD(2) - t97, t4 * qJD(2) + t106, 0, t5 * qJD(2) + t107, t2 * qJD(2) + t100; 0, 0, 0, 0, 0, 0, 0, -t48, -t134, 0, 0, 0, t78, t77, t83, 0, t84, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, -t54, t8 * qJD(2) + t99; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t48, 0, -t54, t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t102; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t18;
