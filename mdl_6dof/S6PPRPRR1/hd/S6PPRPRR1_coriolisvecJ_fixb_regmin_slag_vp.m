% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6PPRPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d5,d6,theta1,theta2,theta4]';
% 
% Output:
% tauc_reg [6x20]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6PPRPRR1_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRPRR1_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:43:38
% EndTime: 2019-03-08 18:43:41
% DurationCPUTime: 1.32s
% Computational Cost: add. (1513->205), mult. (4298->336), div. (0->0), fcn. (3876->14), ass. (0->136)
t72 = sin(pkin(7));
t80 = sin(qJ(3));
t144 = t72 * t80;
t75 = cos(pkin(12));
t76 = cos(pkin(7));
t141 = t75 * t76;
t71 = sin(pkin(12));
t73 = sin(pkin(6));
t83 = cos(qJ(3));
t156 = (t80 * t141 + t71 * t83) * t73;
t77 = cos(pkin(6));
t41 = t77 * t144 + t156;
t62 = t77 * qJD(1) + qJD(2);
t37 = qJD(1) * t156 + t62 * t144;
t70 = sin(pkin(13));
t74 = cos(pkin(13));
t49 = (t70 * t80 - t74 * t83) * t72;
t82 = cos(qJ(5));
t123 = t82 * qJD(3);
t63 = -qJD(6) + t123;
t155 = -qJD(6) - t63;
t78 = sin(qJ(6));
t127 = qJD(6) * t78;
t79 = sin(qJ(5));
t117 = t79 * t127;
t81 = cos(qJ(6));
t124 = t81 * qJD(5);
t68 = t79 ^ 2;
t99 = qJD(3) * t68 - t63 * t82;
t154 = -t63 * t117 - t124 * t99;
t126 = qJD(6) * t81;
t116 = t79 * t126;
t121 = qJD(5) * qJD(6);
t125 = t78 * qJD(5);
t45 = (t82 * t125 + t116) * qJD(3) + t78 * t121;
t34 = t74 * t37;
t132 = qJD(1) * t73;
t118 = t75 * t132;
t145 = t71 * t80;
t143 = t72 * t83;
t55 = t62 * t143;
t36 = t76 * t83 * t118 - t132 * t145 + t55;
t35 = qJD(3) * pkin(3) + t36;
t21 = t70 * t35 + t34;
t19 = qJD(3) * pkin(9) + t21;
t46 = -t72 * t118 + t76 * t62 + qJD(4);
t10 = t82 * t19 + t79 * t46;
t97 = t83 * t141 - t145;
t92 = t97 * t73;
t32 = (qJD(1) * t92 + t55) * qJD(3);
t86 = t37 * qJD(3);
t17 = t74 * t32 - t70 * t86;
t4 = t10 * qJD(5) + t79 * t17;
t153 = t4 * t78;
t152 = t4 * t81;
t151 = t74 * pkin(3);
t114 = t82 * t124;
t131 = qJD(3) * t79;
t115 = t78 * t131;
t44 = qJD(3) * t114 - qJD(6) * t115 + t81 * t121;
t150 = t44 * t78;
t56 = t115 - t124;
t149 = t56 * t63;
t58 = t81 * t131 + t125;
t148 = t58 * t63;
t147 = t63 * t78;
t146 = t63 * t81;
t33 = t70 * t37;
t85 = qJD(3) ^ 2;
t142 = t72 * t85;
t140 = t78 * t82;
t139 = t79 * t56;
t138 = t81 * t82;
t137 = t82 * t45;
t84 = qJD(5) ^ 2;
t136 = t84 * t79;
t135 = t84 * t82;
t22 = t70 * t36 + t34;
t106 = pkin(5) * t79 - pkin(10) * t82;
t60 = t106 * qJD(5);
t134 = t22 - t60;
t133 = -t82 ^ 2 + t68;
t64 = t70 * pkin(3) + pkin(9);
t130 = qJD(5) * t64;
t129 = qJD(5) * t79;
t128 = qJD(5) * t82;
t122 = qJD(3) * qJD(5);
t8 = qJD(5) * pkin(10) + t10;
t113 = t63 * t64 + t8;
t111 = t79 * t122;
t16 = t70 * t32 + t74 * t86;
t20 = t74 * t35 - t33;
t110 = t58 * t129 - t44 * t82;
t18 = -qJD(3) * pkin(4) - t20;
t109 = -qJD(3) * t18 - t17;
t107 = t63 * t116;
t98 = -t82 * pkin(5) - t79 * pkin(10) - pkin(4);
t13 = t98 * qJD(3) - t20;
t1 = t81 * t13 - t78 * t8;
t2 = t78 * t13 + t81 * t8;
t40 = t77 * t143 + t92;
t27 = t70 * t40 + t74 * t41;
t52 = -t73 * t75 * t72 + t77 * t76;
t15 = t82 * t27 + t52 * t79;
t26 = -t74 * t40 + t70 * t41;
t105 = t81 * t15 + t78 * t26;
t104 = -t78 * t15 + t81 * t26;
t103 = t79 * t19 - t82 * t46;
t14 = t79 * t27 - t52 * t82;
t50 = (t70 * t83 + t74 * t80) * t72;
t43 = t82 * t50 + t79 * t76;
t102 = t81 * t43 + t78 * t49;
t101 = -t78 * t43 + t81 * t49;
t42 = t79 * t50 - t82 * t76;
t95 = qJD(3) * t22 - t64 * t84 - t16;
t23 = t74 * t36 - t33;
t94 = qJD(5) * (qJD(3) * (-pkin(4) - t151) + t18 + t23);
t93 = t99 * t78;
t3 = -t103 * qJD(5) + t82 * t17;
t7 = -qJD(5) * pkin(5) + t103;
t88 = qJD(5) * t7 + qJD(6) * t13 - t23 * t63 + t3;
t59 = t106 * qJD(3);
t53 = t98 - t151;
t48 = qJD(3) * t49;
t47 = qJD(3) * t50;
t39 = t41 * qJD(3);
t38 = t40 * qJD(3);
t29 = t43 * qJD(5) - t79 * t48;
t28 = -t42 * qJD(5) - t82 * t48;
t25 = t74 * t38 - t70 * t39;
t24 = t70 * t38 + t74 * t39;
t12 = qJD(3) * t60 + t16;
t11 = t81 * t12;
t6 = t15 * qJD(5) + t79 * t25;
t5 = -t14 * qJD(5) + t82 * t25;
t9 = [0, 0, 0, -t39 * qJD(3), -t38 * qJD(3), t16 * t26 + t17 * t27 - t20 * t24 + t21 * t25, 0, 0, 0, 0, 0, -t6 * qJD(5) + (t26 * t129 - t24 * t82) * qJD(3), -t5 * qJD(5) + (t26 * t128 + t24 * t79) * qJD(3), 0, 0, 0, 0, 0 -(-qJD(6) * t105 + t81 * t24 - t78 * t5) * t63 + t104 * t111 + t6 * t56 + t14 * t45 (qJD(6) * t104 + t78 * t24 + t81 * t5) * t63 - t105 * t111 + t6 * t58 + t14 * t44; 0, 0, 0, -t80 * t142, -t83 * t142, t16 * t49 + t17 * t50 - t20 * t47 - t21 * t48, 0, 0, 0, 0, 0, -t29 * qJD(5) + (t49 * t129 - t47 * t82) * qJD(3), -t28 * qJD(5) + (t49 * t128 + t47 * t79) * qJD(3), 0, 0, 0, 0, 0 -(-qJD(6) * t102 - t78 * t28 + t81 * t47) * t63 + t101 * t111 + t29 * t56 + t42 * t45 (qJD(6) * t101 + t81 * t28 + t78 * t47) * t63 - t102 * t111 + t29 * t58 + t42 * t44; 0, 0, 0, 0 (-t97 * t132 + t36 - t55) * qJD(3), t20 * t22 - t21 * t23 + (-t16 * t74 + t17 * t70) * pkin(3), 0.2e1 * t82 * t111, -0.2e1 * t133 * t122, t135, -t136, 0, t79 * t94 + t95 * t82, -t95 * t79 + t82 * t94, t44 * t81 * t79 + (t114 - t117) * t58 (-t56 * t81 - t58 * t78) * t128 + (-t150 - t45 * t81 + (t56 * t78 - t58 * t81) * qJD(6)) * t79, t110 - t154, t107 + t137 + (-t93 - t139) * qJD(5) (-t63 - t123) * t129 (t53 * t127 + t134 * t81) * t63 + (t113 * t126 + t56 * t130 + t78 * t88 - t11) * t82 + (t7 * t126 - t23 * t56 + t153 + t64 * t45 + (-t64 * t147 + (-t64 * t140 + t81 * t53) * qJD(3) + t1) * qJD(5)) * t79 (t53 * t126 - t134 * t78) * t63 + (t58 * t130 + (-qJD(6) * t113 + t12) * t78 + t88 * t81) * t82 + (-t7 * t127 - t23 * t58 + t152 + t64 * t44 + (-t64 * t146 - (t64 * t138 + t78 * t53) * qJD(3) - t2) * qJD(5)) * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t136, -t135, 0, 0, 0, 0, 0, t107 - t137 + (-t93 + t139) * qJD(5), t110 + t154; 0, 0, 0, 0, 0, 0, -t79 * t85 * t82, t133 * t85, 0, 0, 0, t109 * t79, t109 * t82, -t58 * t146 + t150 (t44 + t149) * t81 + (-t45 + t148) * t78, -t63 * t126 + (t63 * t138 + (-t58 + t125) * t79) * qJD(3), t63 * t127 + (-t63 * t140 + (t56 + t124) * t79) * qJD(3), t63 * t131, -pkin(5) * t45 - t152 + (t103 * t78 + t81 * t59) * t63 - t10 * t56 + (pkin(10) * t146 + t7 * t78) * qJD(6) + (-t1 * t79 + (-pkin(10) * t129 - t7 * t82) * t78) * qJD(3), -pkin(5) * t44 + t153 - (-t103 * t81 + t78 * t59) * t63 - t10 * t58 + (-pkin(10) * t147 + t7 * t81) * qJD(6) + (-t7 * t138 + (-pkin(10) * t124 + t2) * t79) * qJD(3); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t56 ^ 2 + t58 ^ 2, t44 - t149, -t148 - t45, t111, t155 * t2 - t78 * t3 - t7 * t58 + t11, t1 * t155 - t78 * t12 - t81 * t3 + t7 * t56;];
tauc_reg  = t9;
