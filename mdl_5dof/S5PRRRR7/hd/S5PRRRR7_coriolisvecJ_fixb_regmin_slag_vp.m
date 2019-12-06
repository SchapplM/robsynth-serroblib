% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
% 
% Output:
% tauc_reg [5x25]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRR7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:12:59
% EndTime: 2019-12-05 17:13:04
% DurationCPUTime: 1.32s
% Computational Cost: add. (1369->185), mult. (3538->282), div. (0->0), fcn. (2626->8), ass. (0->129)
t96 = sin(qJ(5));
t154 = qJD(5) * t96;
t101 = cos(qJ(4));
t102 = cos(qJ(3));
t99 = sin(qJ(2));
t152 = t99 * qJD(1);
t82 = qJD(2) * pkin(6) + t152;
t130 = pkin(7) * qJD(2) + t82;
t61 = t130 * t102;
t52 = t101 * t61;
t98 = sin(qJ(3));
t60 = t130 * t98;
t53 = qJD(3) * pkin(3) - t60;
t97 = sin(qJ(4));
t121 = -t53 * t97 - t52;
t150 = qJD(2) * t102;
t137 = t101 * t150;
t158 = qJD(2) * t98;
t142 = t97 * t158;
t65 = -t137 + t142;
t170 = t65 * pkin(8);
t16 = -t121 - t170;
t100 = cos(qJ(5));
t58 = t100 * t65;
t67 = -t101 * t158 - t150 * t97;
t33 = t67 * t96 - t58;
t103 = cos(qJ(2));
t147 = t103 * qJD(1);
t91 = -pkin(3) * t102 - pkin(2);
t71 = qJD(2) * t91 - t147;
t42 = t65 * pkin(4) + t71;
t181 = t16 * t154 - t42 * t33;
t156 = qJD(3) * t98;
t45 = -t82 * t156 + (-pkin(7) * t156 + t102 * t147) * qJD(2);
t149 = qJD(3) * t102;
t46 = -t82 * t149 + (-pkin(7) * t149 - t147 * t98) * qJD(2);
t135 = t101 * t46 - t97 * t45;
t109 = qJD(4) * t121 + t135;
t146 = qJD(2) * qJD(3);
t131 = t102 * t146;
t93 = qJD(3) + qJD(4);
t36 = qJD(4) * t137 + t101 * t131 - t142 * t93;
t6 = -t36 * pkin(8) + t109;
t145 = -qJD(4) - qJD(5);
t92 = qJD(3) - t145;
t180 = (-t16 * t92 - t6) * t96 + t181;
t75 = t101 * t98 + t102 * t97;
t179 = t93 * t75;
t37 = t179 * qJD(2);
t122 = t100 * t67 + t65 * t96;
t168 = t122 * t33;
t4 = t122 ^ 2 - t33 ^ 2;
t7 = -qJD(5) * t58 + t100 * t36 + t154 * t67 - t37 * t96;
t1 = -t33 * t92 + t7;
t107 = qJD(5) * t122 - t100 * t37 - t96 * t36;
t2 = -t122 * t92 + t107;
t155 = qJD(4) * t97;
t134 = -t155 * t61 + t97 * t46;
t174 = (qJD(4) * t53 + t45) * t101;
t5 = -t37 * pkin(8) + t134 + t174;
t118 = t100 * t6 + t122 * t42 - t96 * t5;
t178 = -0.2e1 * t146;
t119 = t101 * t102 - t97 * t98;
t63 = t119 * t99;
t104 = qJD(3) ^ 2;
t105 = qJD(2) ^ 2;
t177 = (t104 + t105) * t99;
t111 = t119 * t103;
t148 = qJD(4) * t101;
t172 = pkin(6) + pkin(7);
t138 = qJD(3) * t172;
t76 = t98 * t138;
t77 = t102 * t138;
t78 = t172 * t98;
t79 = t172 * t102;
t176 = qJD(1) * t111 + t101 * t76 + t148 * t78 + t155 * t79 + t77 * t97;
t50 = t97 * t61;
t133 = t101 * t53 - t50;
t64 = t67 * pkin(8);
t15 = t133 + t64;
t112 = t75 * t103;
t120 = -t101 * t79 + t78 * t97;
t175 = qJD(1) * t112 + qJD(4) * t120 - t101 * t77 + t97 * t76;
t173 = qJD(5) - t92;
t114 = pkin(3) * t156 - t152;
t171 = pkin(3) * t92;
t169 = t67 * pkin(4);
t166 = t67 * t65;
t165 = t71 * t67;
t164 = -t101 * t60 - t50;
t143 = pkin(3) * t158;
t72 = qJD(2) * t152 + qJD(3) * t143;
t163 = -t102 ^ 2 + t98 ^ 2;
t162 = qJD(2) * pkin(2);
t161 = t100 * t16;
t160 = t100 * t97;
t159 = t104 * t98;
t157 = qJD(2) * t99;
t153 = t104 * t102;
t140 = -pkin(3) * t93 - t53;
t13 = pkin(4) * t93 + t15;
t139 = -pkin(4) * t92 - t13;
t132 = t60 * t97 - t52;
t129 = pkin(4) * t179 + t114;
t126 = t103 * t178;
t125 = t179 * pkin(8) - qJD(5) * (-pkin(8) * t75 - t101 * t78 - t79 * t97) + t176;
t43 = t93 * t119;
t124 = t43 * pkin(8) + qJD(5) * (pkin(8) * t119 - t120) - t175;
t123 = -t96 * t13 - t161;
t41 = t100 * t75 + t119 * t96;
t40 = -t100 * t119 + t75 * t96;
t117 = qJD(2) * t162;
t116 = t71 * t65 - t134;
t110 = -0.2e1 * qJD(3) * t162;
t90 = pkin(3) * t101 + pkin(4);
t62 = t75 * t99;
t57 = -pkin(4) * t119 + t91;
t48 = t143 - t169;
t24 = -t65 ^ 2 + t67 ^ 2;
t23 = pkin(4) * t37 + t72;
t22 = -t67 * t93 - t37;
t21 = t65 * t93 + t36;
t20 = t64 + t164;
t19 = t132 + t170;
t18 = -qJD(2) * t112 - t63 * t93;
t17 = qJD(2) * t111 - t179 * t99;
t10 = qJD(5) * t41 + t100 * t179 + t96 * t43;
t9 = -qJD(5) * t40 + t100 * t43 - t179 * t96;
t3 = [0, 0, -t105 * t99, -t105 * t103, 0, 0, 0, 0, 0, -t102 * t177 + t126 * t98, t102 * t126 + t177 * t98, 0, 0, 0, 0, 0, -t103 * t37 + t157 * t65 + t18 * t93, -t103 * t36 - t157 * t67 - t17 * t93, 0, 0, 0, 0, 0, (t100 * t18 - t96 * t17 + (-t100 * t63 + t62 * t96) * qJD(5)) * t92 - t33 * t157 + t103 * t107, -(t100 * t17 + t96 * t18 + (-t100 * t62 - t63 * t96) * qJD(5)) * t92 - t122 * t157 - t103 * t7; 0, 0, 0, 0, 0.2e1 * t98 * t131, t163 * t178, t153, -t159, 0, -pkin(6) * t153 + t110 * t98, pkin(6) * t159 + t102 * t110, t36 * t75 - t43 * t67, t119 * t36 + t179 * t67 - t37 * t75 - t43 * t65, t43 * t93, -t179 * t93, 0, t114 * t65 - t119 * t72 + t175 * t93 + t179 * t71 + t91 * t37, -t114 * t67 + t176 * t93 + t91 * t36 + t71 * t43 + t72 * t75, -t122 * t9 + t41 * t7, t10 * t122 + t107 * t41 + t33 * t9 - t40 * t7, t9 * t92, -t10 * t92, 0, t42 * t10 + t23 * t40 - t57 * t107 + (-t100 * t124 + t125 * t96) * t92 - t129 * t33, t23 * t41 + t42 * t9 + t57 * t7 + (t100 * t125 + t124 * t96) * t92 - t129 * t122; 0, 0, 0, 0, -t98 * t105 * t102, t163 * t105, 0, 0, 0, t98 * t117, t102 * t117, -t166, t24, t21, t22, 0, -t132 * t93 - t65 * t143 + t165 + (t140 * t97 - t52) * qJD(4) + t135, t164 * t93 + t67 * t143 + (qJD(4) * t140 - t45) * t101 + t116, t168, t4, t1, t2, 0, -(t100 * t19 - t96 * t20) * t92 + t48 * t33 + (-t101 * t96 - t160) * qJD(4) * t171 + ((-pkin(3) * t160 - t90 * t96) * t92 + t123) * qJD(5) + t118, t48 * t122 + (-t145 * t171 * t97 + t19 * t92 - t6) * t96 + (-qJD(5) * t13 - t5 + (-pkin(3) * t148 - qJD(5) * t90 + t20) * t92) * t100 + t181; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t166, t24, t21, t22, 0, -t121 * t93 + t109 + t165, t133 * t93 + t116 - t174, t168, t4, t1, t2, 0, -(-t15 * t96 - t161) * t92 - t33 * t169 + (t139 * t96 - t161) * qJD(5) + t118, -t122 * t169 + (qJD(5) * t139 + t15 * t92 - t5) * t100 + t180; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t168, t4, t1, t2, 0, t123 * t173 + t118, (-t13 * t173 - t5) * t100 + t180;];
tauc_reg = t3;
