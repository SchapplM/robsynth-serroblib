% Calculate minimal parameter regressor of coriolis matrix for
% S5RRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
% 
% Output:
% cmat_reg [(5*%NQJ)%x20]
%   minimal parameter regressor of coriolis matrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function cmat_reg = S5RRPRP4_coriolismatJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_coriolismatJ_fixb_regmin_slag_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From coriolismat_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:58
% EndTime: 2019-12-31 19:53:01
% DurationCPUTime: 0.90s
% Computational Cost: add. (647->172), mult. (1231->164), div. (0->0), fcn. (833->4), ass. (0->127)
t82 = cos(qJ(4));
t80 = sin(qJ(4));
t48 = -t80 * pkin(4) + t82 * qJ(5);
t98 = -qJ(3) + t48;
t145 = t98 * t82;
t81 = sin(qJ(2));
t150 = t81 * pkin(1);
t35 = -t98 + t150;
t146 = t35 * t82;
t143 = t146 / 0.2e1 - t145 / 0.2e1;
t134 = t80 * qJ(5);
t149 = t82 * pkin(4);
t47 = t134 + t149;
t38 = t47 * t80;
t157 = t143 + t38;
t124 = t80 * qJD(5);
t156 = qJD(4) * t48 + t124;
t78 = t80 ^ 2;
t79 = t82 ^ 2;
t64 = t78 - t79;
t76 = qJD(1) + qJD(2);
t155 = t76 * t64;
t66 = qJ(3) + t150;
t101 = qJ(3) / 0.2e1 + t66 / 0.2e1;
t154 = t101 * t82;
t83 = cos(qJ(2));
t148 = t83 * pkin(1);
t144 = t47 * t82;
t136 = pkin(1) * qJD(2);
t71 = t83 * t136;
t58 = t80 * t71;
t73 = qJD(3) * t80;
t141 = t58 + t73;
t60 = t82 * t71;
t74 = qJD(3) * t82;
t140 = t60 + t74;
t72 = t79 * qJD(5);
t139 = t72 - t74;
t137 = pkin(1) * qJD(1);
t70 = t83 * t137;
t61 = t82 * t70;
t138 = t74 - t61;
t100 = (t78 + t79) * t81;
t102 = -pkin(2) - t148;
t65 = -pkin(7) + t102;
t7 = (t100 * t65 + t35 * t83) * pkin(1);
t135 = t7 * qJD(1);
t10 = t38 + t146;
t132 = t10 * qJD(1);
t28 = t35 * t80;
t11 = -t28 + t144;
t131 = t11 * qJD(1);
t22 = -t102 * t150 - t148 * t66;
t130 = t22 * qJD(1);
t129 = t35 * qJD(1);
t36 = pkin(1) * t100;
t128 = t36 * qJD(1);
t127 = t98 * qJD(2);
t126 = t66 * qJD(1);
t125 = t80 * qJD(4);
t123 = t82 * qJD(1);
t122 = t82 * qJD(4);
t121 = t71 + qJD(3);
t120 = qJ(3) * qJD(2);
t119 = qJ(3) * qJD(4);
t118 = qJD(4) * qJ(5);
t117 = -t60 + t139;
t116 = t81 * t136;
t115 = t81 * t137;
t114 = -t150 / 0.2e1;
t113 = t150 / 0.2e1;
t112 = t47 * t129;
t111 = t35 * t123;
t110 = t80 * t126;
t109 = t66 * t123;
t108 = t65 * t125;
t84 = -pkin(2) - pkin(7);
t107 = t84 * t125;
t106 = t80 * t122;
t105 = t65 * t122;
t104 = t84 * t122;
t103 = t78 / 0.2e1 + t79 / 0.2e1;
t45 = t76 * t82;
t99 = -t124 * t82 + t73;
t14 = t38 - t145;
t56 = t82 * t113;
t4 = t56 - t157;
t97 = qJD(1) * t4 - qJD(2) * t14;
t37 = t98 * t80;
t15 = t37 + t144;
t53 = t80 * t113;
t95 = t144 - t28 / 0.2e1 + t37 / 0.2e1;
t3 = t53 + t95;
t96 = qJD(1) * t3 + qJD(2) * t15;
t12 = (-0.1e1 / 0.2e1 + t103) * t150 + t98;
t94 = qJD(1) * t12 + t127;
t92 = t58 + t99;
t86 = (t134 / 0.2e1 + t149 / 0.2e1) * t150;
t1 = (t98 / 0.2e1 - t35 / 0.2e1) * t47 + t86;
t91 = t1 * qJD(1) + t127 * t47;
t57 = t82 * t114;
t8 = t57 + t143;
t90 = qJD(1) * t8 - t127 * t82;
t54 = t80 * t114;
t17 = t101 * t80 + t54;
t89 = qJD(1) * t17 + t120 * t80;
t18 = t56 - t154;
t88 = qJD(1) * t18 - t120 * t82;
t87 = qJD(4) * t47 - qJD(5) * t82 + qJD(3);
t77 = qJ(3) * qJD(3);
t63 = t76 * qJ(3);
t62 = t66 * qJD(3);
t59 = t80 * t70;
t50 = t64 * qJD(4);
t44 = t76 * t79;
t43 = t76 * t80;
t39 = t76 * t150;
t31 = t80 * t45;
t30 = t36 * qJD(2);
t20 = t56 + t154;
t19 = t54 - (qJ(3) + t66) * t80 / 0.2e1;
t13 = t103 * t150 + t113 - t98;
t9 = t57 - t143;
t6 = t53 - t95;
t5 = t56 + t157;
t2 = t86 + (t35 - t98) * t47 / 0.2e1;
t16 = [0, 0, 0, 0, -t116, -t71, t116, t121, -qJD(2) * t22 + t62, -t106, t50, 0, 0, 0, t122 * t66 + t141, -t125 * t66 + t140, qJD(4) * t10 + t92, -t30, -qJD(4) * t11 + t117, t7 * qJD(2) + t35 * t87; 0, 0, 0, 0, -t39, -t70 - t71, t39, t70 + t121, -t130 + t62 + (-pkin(2) * t81 + qJ(3) * t83) * t136, -t106, t50, 0, 0, 0, qJD(4) * t20 + t141 + t59, qJD(4) * t19 + t140 + t61, qJD(4) * t5 + t59 + t92, -t30 - t128, qJD(4) * t6 + t117 - t61, t135 + t13 * qJD(3) + t2 * qJD(4) + t9 * qJD(5) + (t100 * t84 - t83 * t98) * t136; 0, 0, 0, 0, 0, 0, 0, t76, t76 * t66, 0, 0, 0, 0, 0, t43, t45, t43, 0, -t45, qJD(2) * t13 + t129; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t155, -t125, -t122, 0, qJD(2) * t20 - t108 + t109, qJD(2) * t19 - t105 - t110, qJD(2) * t5 - t108 + t132, -t156, qJD(2) * t6 + t105 - t131, t2 * qJD(2) + t156 * t65 + t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t125, t44, qJD(2) * t9 + t108 - t111; 0, 0, 0, 0, t115, t70, -t115, -t70 + qJD(3), t77 + t130, -t106, t50, 0, 0, 0, -qJD(4) * t18 - t59 + t73, -qJD(4) * t17 + t138, -qJD(4) * t4 - t59 + t99, t128, -qJD(4) * t3 - t138 + t72, -qJD(3) * t12 - qJD(4) * t1 - qJD(5) * t8 - t135; 0, 0, 0, 0, 0, 0, 0, qJD(3), t77, -t106, t50, 0, 0, 0, t119 * t82 + t73, -t119 * t80 + t74, qJD(4) * t14 + t99, 0, -qJD(4) * t15 + t139, -t87 * t98; 0, 0, 0, 0, 0, 0, 0, t76, t63, 0, 0, 0, 0, 0, t43, t45, t43, 0, -t45, -t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, t155, -t125, -t122, 0, -t88 - t107, -t89 - t104, -t97 - t107, -t156, -t96 + t104, t156 * t84 - t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t31, -t125, t44, -t90 + t107; 0, 0, 0, 0, 0, 0, 0, -t76, -t120 - t126, 0, 0, 0, 0, 0, -t43, -t45, -t43, 0, t45, qJD(2) * t12 - t129; 0, 0, 0, 0, 0, 0, 0, -t76, -t63, 0, 0, 0, 0, 0, -t43, -t45, -t43, 0, t45, t94; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t125, -t122, -t125, 0, t122, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t155, 0, 0, 0, qJD(2) * t18 - t109, qJD(2) * t17 + t110, qJD(2) * t4 - t132, 0, qJD(2) * t3 + t131, qJD(2) * t1 - t112; 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, -t155, 0, 0, 0, t88, t89, t97, 0, t96, t91; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(5), qJ(5) * qJD(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJD(4), t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t44, qJD(2) * t8 + t111; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t31, 0, -t44, t90; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJD(4), -t118; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0;];
cmat_reg = t16;
