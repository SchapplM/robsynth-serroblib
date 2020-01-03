% Calculate minimal parameter regressor of coriolis matrix for
% S5RPPRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x23]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:06
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPPRR11_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR11_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:05:56
% EndTime: 2019-12-31 18:05:59
% DurationCPUTime: 0.86s
% Computational Cost: add. (341->137), mult. (831->188), div. (0->0), fcn. (627->4), ass. (0->110)
t62 = sin(qJ(5));
t128 = 0.2e1 * t62;
t63 = sin(qJ(4));
t57 = t63 ^ 2;
t65 = cos(qJ(4));
t59 = t65 ^ 2;
t117 = t57 + t59;
t56 = t62 ^ 2;
t64 = cos(qJ(5));
t58 = t64 ^ 2;
t36 = t58 - t56;
t119 = t64 * t65;
t83 = t119 * t128;
t66 = qJD(1) * t83 - t36 * qJD(4);
t126 = t63 * pkin(7);
t125 = t65 * pkin(4);
t124 = t59 * t64;
t60 = -pkin(6) + qJ(2);
t123 = t60 * t62;
t29 = t125 + t126;
t122 = t62 * t29;
t121 = t64 * t29;
t120 = t64 * t60;
t61 = pkin(1) + qJ(3);
t97 = t65 * qJD(4);
t45 = t62 * t97;
t108 = qJD(5) * t64;
t48 = t63 * t108;
t118 = -t45 - t48;
t35 = t57 - t59;
t78 = t63 * pkin(4) - t65 * pkin(7);
t70 = t78 + t61;
t9 = t123 * t63 - t64 * t70;
t94 = t65 * t123;
t1 = -t9 * t65 + (t94 + t121) * t63;
t116 = t1 * qJD(1);
t95 = t63 * t120;
t10 = t62 * t70 + t95;
t93 = t60 * t119;
t2 = t10 * t65 + (-t93 + t122) * t63;
t115 = t2 * qJD(1);
t5 = -t123 * t59 - t9 * t63;
t114 = t5 * qJD(1);
t6 = -t10 * t63 - t59 * t120;
t113 = t6 * qJD(1);
t112 = qJD(2) * t63;
t111 = qJD(3) * t63;
t110 = qJD(4) * t64;
t109 = qJD(5) * t62;
t81 = 0.1e1 / 0.2e1 + t57 / 0.2e1 + t59 / 0.2e1;
t11 = t81 * t62;
t107 = t11 * qJD(1);
t12 = t81 * t64;
t106 = t12 * qJD(1);
t24 = t117 * t62;
t105 = t24 * qJD(1);
t25 = t35 * t62;
t104 = t25 * qJD(1);
t26 = t117 * t64;
t103 = t26 * qJD(1);
t27 = t35 * t64;
t102 = t27 * qJD(1);
t101 = t35 * qJD(1);
t100 = t61 * qJD(1);
t50 = t63 * qJD(1);
t99 = t63 * qJD(4);
t51 = t65 * qJD(1);
t98 = t65 * qJD(2);
t96 = t65 * qJD(5);
t92 = t62 * t96;
t91 = t64 * t96;
t90 = t62 * t108;
t89 = t62 * t110;
t88 = t62 * t51;
t87 = t64 * t97;
t86 = t64 * t51;
t85 = t63 * t97;
t84 = t63 * t51;
t82 = -qJD(2) + t100;
t79 = qJD(4) * t83;
t47 = t63 * t109;
t77 = t47 - t87;
t76 = (-qJD(5) - t50) * t65;
t75 = t126 / 0.2e1 + t125 / 0.2e1;
t69 = t29 / 0.2e1 + t75;
t8 = t69 * t64;
t74 = pkin(4) * t62 * qJD(4) + t8 * qJD(1);
t7 = t69 * t62;
t73 = pkin(4) * t110 - t7 * qJD(1);
t72 = t64 * t76;
t17 = (t56 / 0.2e1 - t58 / 0.2e1) * t65;
t71 = -t17 * qJD(1) + t89;
t68 = qJD(1) * t124 * t62 + t17 * qJD(4);
t23 = t36 * t59;
t67 = t23 * qJD(1) + t79;
t55 = qJ(2) * qJD(1);
t54 = qJ(2) * qJD(2);
t49 = t97 / 0.2e1;
t46 = t64 * t50;
t44 = t62 * t99;
t43 = t62 * t50;
t22 = t46 + t108;
t21 = t43 + t109;
t20 = (t50 + qJD(5) / 0.2e1) * t65;
t15 = t17 * qJD(5);
t14 = -t124 / 0.2e1 + (-t57 / 0.2e1 + 0.1e1 / 0.2e1) * t64;
t13 = (-0.1e1 / 0.2e1 + t117 / 0.2e1) * t62;
t4 = -t94 + t121 / 0.2e1 - t75 * t64;
t3 = -t93 - t122 / 0.2e1 + t75 * t62;
t16 = [0, 0, 0, 0, qJD(2), t54, qJD(2), qJD(3), t61 * qJD(3) + t54, -t85, t35 * qJD(4), 0, 0, 0, t61 * t97 + t111, qJD(3) * t65 - t61 * t99, -t58 * t85 - t59 * t90, -t23 * qJD(5) + t63 * t79, -t27 * qJD(4) - t63 * t92, t25 * qJD(4) - t63 * t91, t85, -t24 * qJD(2) + t1 * qJD(4) + t6 * qJD(5) + t64 * t111, -t26 * qJD(2) - t2 * qJD(4) - t5 * qJD(5) - t62 * t111; 0, 0, 0, 0, qJD(1), t55, qJD(1), 0, t55, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t105, -t103; 0, 0, 0, 0, 0, 0, 0, qJD(1), t100, 0, 0, 0, 0, 0, t50, t51, 0, 0, 0, 0, 0, t14 * qJD(5) + t46, t13 * qJD(5) - t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, t101, -t99, -t97, 0, t61 * t51 - t60 * t99, -t61 * t50 - t60 * t97, -t15 + (-t58 * t51 - t89) * t63, t63 * t66 - 0.2e1 * t65 * t90, t45 - t102, t87 + t104, t20, t116 + (t78 * t62 - t95) * qJD(4) + t4 * qJD(5), -t115 + (-pkin(7) * t119 + (pkin(4) * t64 + t123) * t63) * qJD(4) + t3 * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t68, -t67, t62 * t76, t72, t49, t14 * qJD(3) + t4 * qJD(4) - t10 * qJD(5) + t113, t13 * qJD(3) + t3 * qJD(4) + t9 * qJD(5) - t114; 0, 0, 0, 0, -qJD(1), -t55, -qJD(1), 0, -qJD(3) - t55, 0, 0, 0, 0, 0, -t97, t99, 0, 0, 0, 0, 0, t77 + t105, t103 - t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t51, t50, 0, 0, 0, 0, 0, -t86, t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t21, t22; 0, 0, 0, 0, 0, 0, 0, -qJD(1), -t82, 0, 0, 0, 0, 0, -t50, -t51, 0, 0, 0, 0, 0, -t12 * qJD(5) - t46, t11 * qJD(5) + t43; 0, 0, 0, 0, 0, 0, 0, 0, qJD(1), 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t99, -t97, 0, 0, 0, 0, 0, -t64 * t99 - t92, t44 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t106 + t118, t77 + t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, -t101, 0, 0, 0, -t82 * t65, t82 * t63, t58 * t84 - t15, t72 * t128, t48 + t102, -t47 - t104, -t20, -t8 * qJD(5) + t64 * t98 - t116, qJD(5) * t7 - t62 * t98 + t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t51, -t50, 0, 0, 0, 0, 0, t86, -t88; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t90, t36 * qJD(5), 0, 0, 0, -pkin(4) * t109, -pkin(4) * t108; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71, -t66, t22, -t21, -t51 / 0.2e1, -pkin(7) * t108 - t74, pkin(7) * t109 - t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t68, t67, (t88 - t110) * t63, t64 * t84 + t44, t49, t12 * qJD(3) + t8 * qJD(4) - t62 * t112 - t113, -qJD(3) * t11 - qJD(4) * t7 - t112 * t64 + t114; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t43, -t46; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t106, -t107; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t71, t66, -t46, t43, t51 / 0.2e1, t74, t73; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
