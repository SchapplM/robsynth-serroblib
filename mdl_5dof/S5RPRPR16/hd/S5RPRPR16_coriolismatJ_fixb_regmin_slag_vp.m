% Calculate minimal parameter regressor of coriolis matrix for
% S5RPRPR16
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x24]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RPRPR16_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPR16_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:39:35
% EndTime: 2019-12-31 18:39:38
% DurationCPUTime: 1.00s
% Computational Cost: add. (480->145), mult. (1018->205), div. (0->0), fcn. (799->4), ass. (0->122)
t75 = sin(qJ(3));
t71 = t75 ^ 2;
t77 = cos(qJ(3));
t73 = t77 ^ 2;
t148 = -t71 / 0.2e1 - t73 / 0.2e1;
t117 = t75 * qJD(4);
t133 = t77 * qJ(4);
t42 = -t75 * pkin(3) + t133;
t89 = t42 * qJD(3) + t117;
t74 = sin(qJ(5));
t70 = t74 ^ 2;
t76 = cos(qJ(5));
t72 = t76 ^ 2;
t45 = t70 - t72;
t139 = t75 * t76;
t98 = 0.2e1 * t74 * t139;
t81 = qJD(1) * t98 - t45 * qJD(3);
t78 = -pkin(3) - pkin(7);
t147 = qJD(3) * (t75 * t78 + t133) + t117;
t79 = -pkin(1) - pkin(6);
t145 = -pkin(4) + t79;
t38 = t145 * t77;
t143 = t38 * t76;
t142 = t71 * t76;
t65 = t75 * qJ(4);
t41 = t77 * pkin(3) + t65;
t29 = t77 * pkin(7) + t41;
t141 = t74 * t29;
t140 = t74 * t75;
t138 = t76 * t29;
t46 = t71 - t73;
t36 = qJ(2) - t42;
t24 = t75 * pkin(7) + t36;
t10 = t74 * t24 + t143;
t1 = -t141 * t77 + (t10 - t143) * t75;
t137 = t1 * qJD(1);
t11 = t76 * t24 - t38 * t74;
t2 = -t11 * t75 + t138 * t77 - t38 * t140;
t136 = t2 * qJD(1);
t37 = t145 * t75;
t3 = -t10 * t77 - t37 * t139;
t135 = t3 * qJD(1);
t4 = -t11 * t77 + t37 * t140;
t134 = t4 * qJD(1);
t132 = qJD(2) * t77;
t131 = qJD(4) * t77;
t130 = qJD(5) * t76;
t129 = qJD(5) * t77;
t128 = qJD(5) * t78;
t12 = t36 * t77 + t41 * t75;
t127 = t12 * qJD(1);
t13 = -t36 * t75 + t41 * t77;
t126 = t13 * qJD(1);
t95 = -0.1e1 / 0.2e1 + t148;
t16 = t95 * t74;
t125 = t16 * qJD(1);
t17 = t95 * t76;
t124 = t17 * qJD(1);
t33 = t46 * t74;
t123 = t33 * qJD(1);
t35 = t46 * t76;
t122 = t35 * qJD(1);
t121 = t36 * qJD(1);
t120 = t46 * qJD(1);
t119 = t73 * qJD(1);
t118 = t74 * qJD(3);
t63 = t75 * qJD(1);
t62 = t75 * qJD(3);
t116 = t75 * qJD(5);
t115 = t76 * qJD(3);
t64 = t77 * qJD(1);
t114 = t77 * qJD(3);
t113 = qJ(2) * qJD(3);
t112 = qJ(4) * qJD(5);
t111 = qJD(1) * qJ(2);
t110 = qJD(3) * qJ(4);
t109 = t74 * t129;
t108 = t76 * t129;
t107 = t41 * t121;
t106 = t36 * t64;
t105 = t75 * t115;
t104 = t74 * t130;
t103 = t74 * t115;
t102 = t75 * t114;
t58 = t75 * t64;
t101 = t79 * t114;
t100 = t75 * t111;
t99 = t77 * t111;
t97 = t76 * t58;
t96 = -qJD(2) + t131;
t93 = qJD(3) * t98;
t91 = -t119 - t129;
t90 = t73 * qJD(4) - t132;
t88 = -t78 * t77 / 0.2e1 + t65 / 0.2e1;
t87 = (-qJD(5) - t64) * t140;
t80 = t29 / 0.2e1 + t88;
t9 = t80 * t76;
t86 = -t9 * qJD(1) + t74 * t110;
t8 = t80 * t74;
t85 = -t8 * qJD(1) - t76 * t110;
t25 = (t72 / 0.2e1 - t70 / 0.2e1) * t75;
t84 = -t25 * qJD(1) + t103;
t83 = t74 * qJD(1) * t142 + t25 * qJD(3);
t34 = t45 * t71;
t82 = t34 * qJD(1) + t93;
t61 = -t62 / 0.2e1;
t60 = t79 * t62;
t59 = t76 * t116;
t57 = t76 * t114;
t56 = t76 * t64;
t55 = t74 * t114;
t54 = t74 * t62;
t53 = t74 * t64;
t32 = -t56 - t130;
t31 = -qJD(5) * t74 - t53;
t30 = t58 + t116 / 0.2e1;
t23 = t25 * qJD(5);
t19 = t142 / 0.2e1 + (t73 / 0.2e1 - 0.1e1 / 0.2e1) * t76;
t18 = (-0.1e1 / 0.2e1 - t148) * t74;
t6 = -t37 * t74 - t138 / 0.2e1 + t88 * t76;
t5 = t37 * t76 - t141 / 0.2e1 + t88 * t74;
t7 = [0, 0, 0, 0, qJD(2), qJ(2) * qJD(2), -t102, t46 * qJD(3), 0, 0, 0, qJD(2) * t75 + t113 * t77, -t113 * t75 + t132, 0, -t12 * qJD(3) + t75 * t96, -t13 * qJD(3) + t90, (qJD(3) * t41 - t96) * t36, t102 * t70 + t104 * t71, -t34 * qJD(5) + t77 * t93, -t33 * qJD(3) + t108 * t75, -t35 * qJD(3) - t109 * t75, -t102, t1 * qJD(3) + t4 * qJD(5) + t74 * t90, -t2 * qJD(3) - t3 * qJD(5) + t76 * t90; 0, 0, 0, 0, qJD(1), t111, 0, 0, 0, 0, 0, t63, t64, 0, -t63, -t64, t121, 0, 0, 0, 0, 0, t18 * qJD(5) - t53, t19 * qJD(5) - t56; 0, 0, 0, 0, 0, 0, -t58, t120, -t62, -t114, 0, -t60 + t99, -t100 - t101, -t89, t60 - t127, t101 - t126, t89 * t79 + t107, t23 + (t63 * t70 + t103) * t77, -0.2e1 * t75 * t104 + t81 * t77, -t105 - t123, t54 - t122, -t30, t5 * qJD(5) + t38 * t118 - t147 * t76 + t137, t6 * qJD(5) + t38 * t115 + t147 * t74 - t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, t58, t119, t60 - t106, 0, 0, 0, 0, 0, t119 * t74 - t105, t119 * t76 + t54; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83, -t82, t59 + t97, t87, t61, t18 * qJD(2) + t5 * qJD(3) - t11 * qJD(5) + t134, t19 * qJD(2) + t6 * qJD(3) + t10 * qJD(5) - t135; 0, 0, 0, 0, -qJD(1), -t111, 0, 0, 0, 0, 0, -t63, -t64, 0, t63, t64, -t121, 0, 0, 0, 0, 0, -t16 * qJD(5) + t53, -t17 * qJD(5) + t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t62, -t114, 0, t62, t114, t89, 0, 0, 0, 0, 0, t59 + t55, -t116 * t74 + t57; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t62, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t105 + t109 - t125, t108 - t54 - t124; 0, 0, 0, 0, 0, 0, t58, -t120, 0, 0, 0, -t99, t100, 0, t127, t126, -t107, -t58 * t70 + t23, 0.2e1 * t76 * t87, -t109 + t123, -t108 + t122, t30, t8 * qJD(5) - t137, t9 * qJD(5) + t136; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), qJ(4) * qJD(4), -t104, t45 * qJD(5), 0, 0, 0, qJD(4) * t74 + t112 * t76, qJD(4) * t76 - t112 * t74; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(3), t110, 0, 0, 0, 0, 0, t118, t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t84, -t81, t31, t32, t63 / 0.2e1, -t74 * t128 - t85, -t76 * t128 - t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t58, -t119, t106, 0, 0, 0, 0, 0, t91 * t74, t91 * t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(3), -t110, 0, 0, 0, 0, 0, -t118, -t115; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t82, t55 - t97, t58 * t74 + t57, t61, t16 * qJD(2) - t8 * qJD(3) + t74 * t131 - t134, t17 * qJD(2) - t9 * qJD(3) + t76 * t131 + t135; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, t124; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t84, t81, t53, t56, -t63 / 0.2e1, t85, t86; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t53, t56; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t7;
