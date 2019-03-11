% Calculate minimal parameter regressor of joint inertia matrix time derivative for
% S6RPRPRR12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
% 
% Output:
% MMD_reg [((6+1)*6/2)x31]
%   minimal parameter regressor of inerta matrix time derivative
%   (only lower left triangular matrix (including diagonal) due to symmetry

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MMD_reg = S6RPRPRR12_inertiaDJ_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_inertiaDJ_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_inertiaDJ_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_inertiaDJ_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From inertiaD_joint_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:19:57
% EndTime: 2019-03-09 04:20:01
% DurationCPUTime: 1.38s
% Computational Cost: add. (908->170), mult. (2014->285), div. (0->0), fcn. (1668->6), ass. (0->115)
t109 = qJD(5) + qJD(6);
t70 = sin(qJ(6));
t71 = sin(qJ(5));
t73 = cos(qJ(6));
t74 = cos(qJ(5));
t42 = t70 * t74 + t73 * t71;
t140 = t109 * t42;
t72 = sin(qJ(3));
t142 = t140 * t72;
t64 = t72 * qJD(3);
t75 = cos(qJ(3));
t126 = t73 * t74;
t83 = t70 * t71 - t126;
t13 = t140 * t75 - t83 * t64;
t25 = t42 * t72;
t24 = t83 * t72;
t122 = t75 * qJ(4);
t139 = -t72 * pkin(3) + t122;
t44 = qJ(2) - t139;
t35 = t72 * pkin(8) + t44;
t77 = -pkin(1) - pkin(7);
t48 = (pkin(4) - t77) * t75;
t36 = t71 * t48;
t125 = t74 * t35 + t36;
t118 = qJD(5) * t75;
t103 = t74 * t118;
t138 = t71 * t64 - t103;
t67 = t72 ^ 2;
t69 = t75 ^ 2;
t91 = (t67 - t69) * qJD(3);
t66 = t71 ^ 2;
t124 = -t74 ^ 2 + t66;
t90 = t124 * qJD(5);
t113 = t72 * qJD(4);
t62 = t72 * t77;
t46 = -t72 * pkin(4) + t62;
t114 = t46 * qJD(5);
t76 = -pkin(3) - pkin(8);
t137 = (t72 * t76 + t122) * qJD(3) + t113 - t114;
t136 = 2 * qJD(2);
t135 = 0.2e1 * qJD(4);
t134 = pkin(9) * t72;
t120 = qJD(5) * t71;
t112 = t75 * qJD(3);
t99 = t74 * t112;
t30 = -t72 * t120 + t99;
t119 = qJD(5) * t74;
t96 = pkin(3) * t112 + qJ(4) * t64 + qJD(2);
t21 = (qJD(3) * pkin(8) - qJD(4)) * t75 + t96;
t58 = t77 * t64;
t38 = -pkin(4) * t64 + t58;
t6 = -t48 * t119 + t35 * t120 - t74 * t21 - t71 * t38;
t5 = t30 * pkin(9) - t6;
t133 = t73 * t5;
t132 = t75 * pkin(5);
t131 = pkin(9) - t76;
t16 = t74 * t134 + t125;
t130 = t70 * t16;
t128 = t71 * t75;
t127 = t73 * t16;
t121 = qJD(3) * t46;
t117 = qJD(5) * t76;
t116 = qJD(6) * t70;
t115 = qJD(6) * t73;
t111 = qJ(2) * qJD(3);
t110 = qJ(4) * qJD(5);
t108 = pkin(5) * t64;
t107 = pkin(5) * t116;
t106 = pkin(5) * t115;
t105 = t42 * t112;
t104 = t71 * t118;
t101 = t71 * t112;
t100 = t74 * t64;
t98 = t71 * t119;
t97 = t72 * t112;
t92 = -t71 * t21 + t74 * t38;
t93 = -t35 - t134;
t4 = (-pkin(5) * t72 - pkin(9) * t128) * qJD(3) + (t93 * t74 - t36) * qJD(5) + t92;
t95 = t73 * t4 - t70 * t5;
t47 = t131 * t74;
t37 = t74 * t48;
t15 = t93 * t71 + t132 + t37;
t94 = -t15 - t132;
t89 = qJD(5) * (t67 + t69);
t88 = t71 * t99;
t59 = t77 * t112;
t39 = -pkin(4) * t112 + t59;
t87 = t73 * t15 - t130;
t86 = t70 * t15 + t127;
t45 = t131 * t71;
t85 = -t73 * t45 - t70 * t47;
t84 = -t70 * t45 + t73 * t47;
t18 = t109 * t126 - t71 * t116 - t70 * t120;
t80 = -t18 * t75 + t42 * t64;
t79 = t39 + (qJ(4) * t72 - t75 * t76) * qJD(5);
t27 = t139 * qJD(3) + t113;
t60 = t71 * pkin(5) + qJ(4);
t54 = pkin(5) * t119 + qJD(4);
t53 = -0.2e1 * t97;
t41 = qJD(5) * t47;
t40 = t131 * t120;
t32 = t72 * t119 + t101;
t31 = t100 + t104;
t26 = t62 + (-pkin(5) * t74 - pkin(4)) * t72;
t22 = -t75 * qJD(4) + t96;
t19 = -t30 * pkin(5) + t39;
t12 = t70 * t101 - t73 * t99 + t142;
t11 = -t74 * t75 * t115 + (t109 * t128 + t100) * t70 + t138 * t73;
t10 = -t109 * t24 + t105;
t9 = -t85 * qJD(6) + t73 * t40 + t70 * t41;
t8 = t84 * qJD(6) - t70 * t40 + t73 * t41;
t7 = -t125 * qJD(5) + t92;
t2 = -t86 * qJD(6) + t95;
t1 = -t87 * qJD(6) - t70 * t4 - t133;
t3 = [0, 0, 0, 0, t136, qJ(2) * t136, t53, 0.2e1 * t91, 0, 0, 0, 0.2e1 * qJD(2) * t72 + 0.2e1 * t75 * t111, 0.2e1 * qJD(2) * t75 - 0.2e1 * t72 * t111, 0, -0.2e1 * t44 * t112 - 0.2e1 * t22 * t72, -0.2e1 * t22 * t75 + 0.2e1 * t44 * t64, 0.2e1 * t44 * t22, 0.2e1 * t66 * t97 + 0.2e1 * t67 * t98, -0.2e1 * t67 * t90 + 0.4e1 * t72 * t88, 0.2e1 * t72 * t103 - 0.2e1 * t71 * t91, -0.2e1 * t72 * t104 - 0.2e1 * t74 * t91, t53, 0.2e1 * (-t74 * t121 + t7) * t75 + 0.2e1 * (-(-t71 * t35 + t37) * qJD(3) - t39 * t74 + t71 * t114) * t72, 0.2e1 * (t71 * t121 + t6) * t75 + 0.2e1 * (t125 * qJD(3) + t74 * t114 + t39 * t71) * t72, 0.2e1 * t25 * t10, -0.2e1 * t10 * t24 - 0.2e1 * t25 * t12, 0.2e1 * t10 * t75 - 0.2e1 * t25 * t64, -0.2e1 * t12 * t75 + 0.2e1 * t24 * t64, t53, 0.2e1 * t26 * t12 + 0.2e1 * t19 * t24 + 0.2e1 * t2 * t75 - 0.2e1 * t87 * t64, 0.2e1 * t1 * t75 + 0.2e1 * t26 * t10 + 0.2e1 * t19 * t25 + 0.2e1 * t86 * t64; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t71 * t89, t74 * t89, 0, 0, 0, 0, 0, t72 * t12 + t13 * t75, t72 * t10 - t11 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t112, 0, -t58, -t59, -t27, t58, t59, t27 * t77, -t72 * t90 + t88, -t124 * t112 - 0.4e1 * t72 * t98, -t31, t138, 0, -t137 * t74 + t79 * t71, t137 * t71 + t79 * t74, -t10 * t83 - t140 * t25, -t10 * t42 + t12 * t83 + t140 * t24 - t25 * t18, -t13, t80, 0, t60 * t12 + t26 * t18 + t19 * t42 + t54 * t24 + t84 * t64 + t9 * t75, t60 * t10 - t140 * t26 - t19 * t83 + t54 * t25 + t85 * t64 + t8 * t75; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, -t112, 0, t64, t112, t27, 0, 0, 0, 0, 0, t32, t30, 0, 0, 0, 0, 0, t72 * t18 + t105, -t112 * t83 - t142; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135, qJ(4) * t135, -0.2e1 * t98, 0.2e1 * t90, 0, 0, 0, 0.2e1 * qJD(4) * t71 + 0.2e1 * t74 * t110, 0.2e1 * qJD(4) * t74 - 0.2e1 * t71 * t110, 0.2e1 * t83 * t140, 0.2e1 * t140 * t42 + 0.2e1 * t18 * t83, 0, 0, 0, 0.2e1 * t60 * t18 + 0.2e1 * t54 * t42, -0.2e1 * t140 * t60 - 0.2e1 * t54 * t83; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t64, 0, 0, t58, 0, 0, 0, 0, 0, -t31, t138, 0, 0, 0, 0, 0, -t13, t80; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t64, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t32, t30, -t64, t7, t6, 0, 0, t10, -t12, -t64, -t73 * t108 + (t94 * t70 - t127) * qJD(6) + t95, -t133 + (-t4 + t108) * t70 + (t94 * t73 + t130) * qJD(6); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t138, 0, 0, 0, 0, 0, t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, -t71 * t117, -t74 * t117, 0, 0, -t140, -t18, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t120, -t119, 0, 0, 0, 0, 0, -t140, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -0.2e1 * t107, -0.2e1 * t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t10, -t12, -t64, t2, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t13, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t18, 0, t9, t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t140, -t18; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t107, -t106; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
MMD_reg  = t3;
