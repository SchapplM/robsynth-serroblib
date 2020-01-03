% Calculate inertial parameters regressor of coriolis matrix for
% S5RRPPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
% 
% Output:
% cmat_reg [(5*5)x(5*10)]
%   inertial parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPPR4_coriolismatJ_fixb_reg2_slag_vp(qJ, qJD, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_reg2_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_reg2_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR4_coriolismatJ_fixb_reg2_slag_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:27:53
% EndTime: 2019-12-31 19:27:55
% DurationCPUTime: 1.25s
% Computational Cost: add. (965->144), mult. (1687->172), div. (0->0), fcn. (1333->6), ass. (0->112)
t90 = sin(pkin(8));
t91 = cos(pkin(8));
t93 = sin(qJ(2));
t95 = cos(qJ(2));
t58 = (t90 * t93 + t91 * t95) * pkin(1);
t141 = -t58 / 0.2e1;
t137 = t95 * pkin(1);
t107 = -pkin(2) - t137;
t76 = -pkin(3) + t107;
t138 = t93 * pkin(1);
t77 = qJ(3) + t138;
t34 = t76 * t91 - t77 * t90;
t31 = pkin(4) - t34;
t96 = -pkin(2) - pkin(3);
t63 = -qJ(3) * t90 + t91 * t96;
t60 = pkin(4) - t63;
t111 = t60 / 0.2e1 + t31 / 0.2e1;
t151 = t141 - t111;
t86 = qJD(1) + qJD(2);
t94 = cos(qJ(5));
t145 = t86 * t94;
t54 = t90 * t145;
t92 = sin(qJ(5));
t146 = t86 * t92;
t53 = t90 * t146;
t150 = t91 * t145;
t149 = t91 * t146;
t88 = t92 ^ 2;
t89 = t94 ^ 2;
t132 = t88 + t89;
t59 = t132 * t91;
t148 = t86 * t59;
t73 = t89 - t88;
t147 = t86 * t73;
t57 = t137 * t90 - t138 * t91;
t119 = t57 * qJD(1);
t82 = t90 * qJD(3);
t144 = t82 - t119;
t142 = t57 / 0.2e1;
t15 = t132 * t58;
t56 = t59 * qJD(3);
t136 = -qJD(2) * t15 - t56;
t49 = t57 * qJD(2);
t71 = t94 * t82;
t135 = t49 * t94 + t71;
t134 = t49 + t82;
t83 = t91 * qJD(3);
t133 = qJD(2) * t58 + t83;
t35 = t76 * t90 + t77 * t91;
t64 = t91 * qJ(3) + t90 * t96;
t131 = pkin(1) * qJD(1);
t130 = pkin(1) * qJD(2);
t32 = -pkin(7) + t35;
t5 = t15 * t32 + t31 * t57;
t129 = t5 * qJD(1);
t6 = t31 * t90 + t32 * t59;
t128 = t6 * qJD(1);
t7 = -t34 * t57 + t35 * t58;
t127 = t7 * qJD(1);
t102 = t141 + t111;
t8 = t102 * t92;
t126 = t8 * qJD(1);
t9 = t102 * t94;
t125 = t9 * qJD(1);
t124 = qJD(1) * t31;
t123 = qJD(2) * t60;
t13 = -t34 * t90 + t35 * t91;
t122 = t13 * qJD(1);
t121 = t15 * qJD(1);
t33 = (t107 * t93 + t77 * t95) * pkin(1);
t120 = t33 * qJD(1);
t118 = t58 * qJD(1);
t85 = t92 * qJD(5);
t117 = t94 * qJD(5);
t81 = t95 * t130;
t116 = t81 + qJD(3);
t115 = t93 * t130;
t114 = t94 * t119;
t113 = -t57 * t91 / 0.2e1;
t112 = t58 * t90 / 0.2e1;
t110 = -t63 / 0.2e1 - t34 / 0.2e1;
t109 = t64 / 0.2e1 + t35 / 0.2e1;
t108 = t89 / 0.2e1 + t88 / 0.2e1;
t61 = -pkin(7) + t64;
t12 = t59 * t61 + t90 * t60;
t97 = (t32 + t61) * t108;
t2 = (-t108 * t58 + t111) * t90 + (t142 + t97) * t91;
t101 = -qJD(1) * t2 - qJD(2) * t12;
t16 = -t90 * t63 + t91 * t64;
t4 = (t142 + t109) * t91 + (t141 + t110) * t90;
t100 = qJD(1) * t4 + qJD(2) * t16;
t99 = -t83 + t124;
t98 = -t83 + t123;
t87 = qJ(3) * qJD(3);
t80 = t95 * t131;
t79 = t93 * t131;
t74 = t92 * t117;
t72 = t86 * qJ(3);
t70 = t77 * qJD(3);
t69 = t73 * qJD(5);
t66 = t86 * t91;
t65 = t86 * t90;
t62 = -t79 - t115;
t55 = t92 * t145;
t28 = t85 * t91 - t54;
t27 = t117 * t91 + t53;
t21 = (-0.1e1 + t132) * t90 * t83;
t11 = t94 * t151;
t10 = t92 * t151;
t3 = t109 * t91 + t110 * t90 + t112 + t113;
t1 = t111 * t90 + t112 * t132 + t97 * t91 + t113;
t14 = [0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t115, -t81, 0, 0, 0, 0, 0, 0, 0, 0, -t115, 0, t116, qJD(2) * t33 + t70, 0, 0, 0, 0, 0, 0, t134, t133, 0, qJD(2) * t7 + qJD(3) * t13, t74, t69, 0, -t74, 0, 0, -t31 * t85 + t135, -t117 * t31 - t134 * t92, t136, qJD(2) * t5 + qJD(3) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, -t80 - t81, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, t80 + t116, t120 + t70 + (-pkin(2) * t93 + qJ(3) * t95) * t130, 0, 0, 0, 0, 0, 0, t119 + t134, t118 + t133, 0, t127 + (-t57 * t63 + t58 * t64) * qJD(2) + t3 * qJD(3), t74, t69, 0, -t74, 0, 0, qJD(5) * t10 + t114 + t135, t11 * qJD(5) + (-t57 * t86 - t82) * t92, -t121 + t136, t129 + (t15 * t61 + t57 * t60) * qJD(2) + t1 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t86 * t77, 0, 0, 0, 0, 0, 0, t65, t66, 0, qJD(2) * t3 + t122, 0, 0, 0, 0, 0, 0, t54, -t53, -t148, qJD(2) * t1 + t128 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t147, -t117, -t55, t85, 0, qJD(2) * t10 - t117 * t32 - t124 * t92, qJD(2) * t11 - t124 * t94 + t32 * t85, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t79, t80, 0, 0, 0, 0, 0, 0, 0, 0, t79, 0, -t80 + qJD(3), t87 - t120, 0, 0, 0, 0, 0, 0, t144, t83 - t118, 0, qJD(3) * t4 - t127, t74, t69, 0, -t74, 0, 0, -qJD(5) * t8 - t114 + t71, -t9 * qJD(5) - t144 * t92, -t56 + t121, qJD(3) * t2 - t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t87, 0, 0, 0, 0, 0, 0, t82, t83, 0, t16 * qJD(3), t74, t69, 0, -t74, 0, 0, -t60 * t85 + t71, -t117 * t60 - t82 * t92, -t56, t12 * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t86, t72, 0, 0, 0, 0, 0, 0, t65, t66, 0, t100, 0, 0, 0, 0, 0, 0, t54, -t53, -t148, -t101 + t21; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t55, t147, -t117, -t55, t85, 0, -t117 * t61 - t123 * t92 - t126, -t123 * t94 + t61 * t85 - t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -qJ(3) * qJD(2) - qJD(1) * t77, 0, 0, 0, 0, 0, 0, -t65, -t66, 0, -qJD(2) * t4 - t122, 0, 0, 0, 0, 0, 0, t28, t27, t148, -qJD(2) * t2 - t128; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t86, -t72, 0, 0, 0, 0, 0, 0, -t65, -t66, 0, -t100, 0, 0, 0, 0, 0, 0, t28, t27, t148, t101; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t117 * t90 + t149, t85 * t90 + t150, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t85, -t117, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t147, 0, t55, 0, 0, t8 * qJD(2) + t92 * t99, t9 * qJD(2) + t94 * t99, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55, -t147, 0, t55, 0, 0, t92 * t98 + t126, t94 * t98 + t125, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, -t150, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t14;
