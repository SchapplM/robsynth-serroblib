% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta4]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x26]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRP8_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP8_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP8_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP8_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:26:04
% EndTime: 2019-03-09 03:26:07
% DurationCPUTime: 1.04s
% Computational Cost: add. (1694->154), mult. (3327->272), div. (0->0), fcn. (3038->6), ass. (0->99)
t58 = sin(qJ(5));
t111 = qJD(5) * t58;
t59 = sin(qJ(3));
t61 = cos(qJ(3));
t112 = sin(pkin(9));
t89 = qJD(3) * t112;
t113 = cos(pkin(9));
t90 = qJD(3) * t113;
t34 = -t59 * t90 - t61 * t89;
t93 = t113 * t61;
t39 = -t112 * t59 + t93;
t125 = t39 * t34;
t33 = t59 * t89 - t61 * t90;
t92 = t112 * t61;
t40 = -t113 * t59 - t92;
t127 = t33 * t40;
t138 = 0.2e1 * t125 + 0.2e1 * t127;
t37 = t39 ^ 2;
t38 = t40 ^ 2;
t60 = cos(qJ(5));
t139 = (t38 + t37) * t111 - t138 * t60;
t137 = (-t112 * t33 + t113 * t34) * pkin(3);
t115 = t59 * pkin(3) + qJ(2);
t22 = -t40 * pkin(4) - t39 * pkin(8) + t115;
t62 = -pkin(1) - pkin(7);
t114 = qJ(4) - t62;
t42 = t114 * t59;
t28 = -t113 * t42 - t114 * t92;
t136 = t58 * t22 + t60 * t28;
t106 = t61 * qJD(3);
t31 = -t59 * qJD(4) - t114 * t106;
t108 = t59 * qJD(3);
t69 = -t61 * qJD(4) + t114 * t108;
t13 = t112 * t69 + t113 * t31;
t43 = pkin(3) * t106 + qJD(2);
t14 = -t33 * pkin(4) - t34 * pkin(8) + t43;
t4 = -qJD(5) * t136 - t58 * t13 + t60 * t14;
t110 = t40 * qJD(6);
t116 = t33 * qJ(6);
t54 = qJD(5) * t60;
t3 = t28 * t111 - t60 * t13 - t58 * t14 - t22 * t54;
t1 = -t110 - t3 - t116;
t129 = pkin(5) * t33;
t2 = t129 - t4;
t6 = -t40 * qJ(6) + t136;
t82 = t60 * t22 - t58 * t28;
t7 = t40 * pkin(5) - t82;
t85 = t58 * t6 - t60 * t7;
t134 = t85 * qJD(5) - t1 * t60 - t2 * t58;
t86 = t58 * t7 + t6 * t60;
t133 = qJD(5) * t86 + t1 * t58 - t2 * t60;
t84 = t60 * pkin(5) + t58 * qJ(6);
t65 = t84 * qJD(5) - t60 * qJD(6);
t132 = 0.2e1 * qJD(2);
t131 = 0.2e1 * qJD(5);
t130 = 0.2e1 * qJD(6);
t12 = t112 * t31 - t113 * t69;
t128 = t12 * t58;
t50 = t112 * pkin(3) + pkin(8);
t126 = t33 * t50;
t124 = t39 * t60;
t123 = t40 * t50;
t122 = t58 * t33;
t121 = t58 * t34;
t120 = t60 * t34;
t56 = t58 ^ 2;
t57 = t60 ^ 2;
t118 = t56 - t57;
t117 = t56 + t57;
t109 = t58 * qJD(6);
t105 = qJ(2) * qJD(3);
t104 = -0.2e1 * t122;
t51 = -t113 * pkin(3) - pkin(4);
t102 = t51 * t131;
t101 = -0.2e1 * t39 * t121 - t37 * t54;
t100 = t50 * t111;
t99 = t40 * t54;
t98 = t50 * t54;
t97 = t58 * t54;
t96 = t117 * t33;
t95 = -0.4e1 * t58 * t124;
t91 = t118 * qJD(5);
t83 = pkin(5) * t58 - qJ(6) * t60;
t32 = -pkin(5) * t111 + qJ(6) * t54 + t109;
t35 = -t84 + t51;
t80 = t39 * t32 - t34 * t35;
t78 = t34 * t51 + t126;
t77 = t39 * t51 + t123;
t27 = -t112 * t42 + t114 * t93;
t76 = t99 + t122;
t75 = t40 * t111 - t60 * t33;
t74 = t39 * t54 + t121;
t73 = -t39 * t111 + t120;
t5 = t83 * t34 + t65 * t39 + t12;
t70 = -t5 + (t35 * t39 + t123) * qJD(5);
t67 = t12 * t39 + t13 * t40 + t27 * t34 + t28 * t33;
t8 = t83 * t39 + t27;
t66 = qJD(5) * t8 + t126 - t80;
t9 = [0, 0, 0, 0, t132, qJ(2) * t132, -0.2e1 * t59 * t106, 0.2e1 * (t59 ^ 2 - t61 ^ 2) * qJD(3), 0, 0, 0, 0.2e1 * qJD(2) * t59 + 0.2e1 * t61 * t105, 0.2e1 * qJD(2) * t61 - 0.2e1 * t59 * t105, 0.2e1 * t67, 0.2e1 * t115 * t43 + 0.2e1 * t27 * t12 + 0.2e1 * t28 * t13, 0.2e1 * t57 * t125 - 0.2e1 * t37 * t97, t118 * t37 * t131 + t34 * t95, -0.2e1 * t40 * t120 + 0.2e1 * t75 * t39, 0.2e1 * t40 * t121 + 0.2e1 * t76 * t39, 0.2e1 * t127, 0.2e1 * t39 * t128 + 0.2e1 * t74 * t27 - 0.2e1 * t82 * t33 - 0.2e1 * t4 * t40, 0.2e1 * t12 * t124 + 0.2e1 * t136 * t33 + 0.2e1 * t73 * t27 - 0.2e1 * t3 * t40, 0.2e1 * t8 * t121 + 0.2e1 * t2 * t40 + 0.2e1 * t7 * t33 + 0.2e1 * (t5 * t58 + t8 * t54) * t39, -0.2e1 * t133 * t39 - 0.2e1 * t85 * t34, -0.2e1 * t8 * t120 - 0.2e1 * t1 * t40 - 0.2e1 * t6 * t33 + 0.2e1 * (t111 * t8 - t5 * t60) * t39, 0.2e1 * t6 * t1 + 0.2e1 * t7 * t2 + 0.2e1 * t8 * t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t138, -t67, 0, 0, 0, 0, 0, t40 * t104 - t38 * t54 + t101, t139 (-t99 + t104) * t40 + t101, 0, -t139, t134 * t40 - t86 * t33 - t8 * t34 - t5 * t39; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t117 * t127 + 0.2e1 * t125; 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t106, 0, -t62 * t108, -t62 * t106, -t137 (t112 * t13 - t113 * t12) * pkin(3), t58 * t120 - t39 * t91, qJD(5) * t95 - t118 * t34, -t76, t75, 0, -t12 * t60 + t78 * t58 + (t27 * t58 + t77 * t60) * qJD(5), t128 + t78 * t60 + (t27 * t60 - t77 * t58) * qJD(5), t66 * t58 + t70 * t60, -t134, t58 * t70 - t60 * t66, -t134 * t50 - t8 * t32 + t5 * t35; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t108, -t106, 0, t137, 0, 0, 0, 0, 0, t73, -t74, t73, -t96, t74, -t50 * t96 + t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0.2e1 * t97, -0.2e1 * t91, 0, 0, 0, t58 * t102, t60 * t102, 0.2e1 * t35 * t111 + 0.2e1 * t32 * t60, 0, 0.2e1 * t32 * t58 - 0.2e1 * t35 * t54, -0.2e1 * t35 * t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t43, 0, 0, 0, 0, 0, t75, t76, t75, -t117 * t34, -t76, t133; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t73, -t74, -t33, t4, t3, t4 - 0.2e1 * t129, -t84 * t34 + (qJD(5) * t83 - t109) * t39, -0.2e1 * t110 - t3 - 0.2e1 * t116, -t2 * pkin(5) + t1 * qJ(6) + t6 * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t76, -t75, t76, 0, t75, t33 * t83 + t40 * t65; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, -t111, 0, -t98, t100, -t98, -t65, -t100, -t65 * t50; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t111, -t54, -t111, 0, t54, t32; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, qJ(6) * t130; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t33, t73, 0, t2; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, 0, t98; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t9;
