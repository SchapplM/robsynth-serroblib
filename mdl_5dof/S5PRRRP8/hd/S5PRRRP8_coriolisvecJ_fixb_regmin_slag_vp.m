% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
% 
% Output:
% tauc_reg [5x22]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5PRRRP8_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:00:36
% EndTime: 2019-12-05 17:00:42
% DurationCPUTime: 1.55s
% Computational Cost: add. (1582->266), mult. (4041->376), div. (0->0), fcn. (2794->8), ass. (0->140)
t85 = cos(qJ(3));
t132 = t85 * qJD(2);
t72 = -qJD(4) + t132;
t131 = qJD(2) * qJD(3);
t82 = sin(qJ(3));
t112 = t82 * t131;
t105 = pkin(4) * t112;
t84 = cos(qJ(4));
t136 = qJD(4) * t84;
t81 = sin(qJ(4));
t137 = qJD(4) * t81;
t79 = sin(pkin(5));
t141 = qJD(2) * t79;
t86 = cos(qJ(2));
t123 = t86 * t141;
t139 = qJD(3) * t82;
t142 = qJD(1) * t85;
t143 = qJD(1) * t79;
t83 = sin(qJ(2));
t118 = t83 * t143;
t65 = qJD(2) * pkin(7) + t118;
t80 = cos(pkin(5));
t22 = -t65 * t139 + (qJD(3) * t80 + t123) * t142;
t158 = t80 * t82;
t71 = qJD(1) * t158;
t45 = t85 * t65 + t71;
t36 = qJD(3) * pkin(8) + t45;
t103 = pkin(3) * t82 - pkin(8) * t85;
t63 = t103 * qJD(3);
t43 = (t63 + t118) * qJD(2);
t117 = t86 * t143;
t68 = -t85 * pkin(3) - t82 * pkin(8) - pkin(2);
t46 = qJD(2) * t68 - t117;
t111 = -t36 * t136 - t46 * t137 - t81 * t22 + t84 * t43;
t2 = -t105 - t111;
t10 = t84 * t36 + t81 * t46;
t7 = -t72 * qJ(5) + t10;
t178 = t7 * t72 + t2;
t154 = t85 * t86;
t177 = -(t81 * t154 - t83 * t84) * t143 + t68 * t137 - t84 * t63;
t176 = -(t84 * t154 + t81 * t83) * t143 + t68 * t136 + t81 * t63;
t44 = t80 * t142 - t82 * t65;
t119 = t82 * t136;
t130 = qJD(3) * qJD(4);
t134 = t81 * qJD(3);
t40 = qJD(2) * (t85 * t134 + t119) + t81 * t130;
t140 = qJD(2) * t82;
t60 = t84 * t140 + t134;
t175 = t60 ^ 2;
t106 = t82 * t117;
t138 = qJD(3) * t85;
t23 = qJD(2) * t106 + qJD(3) * t71 + t65 * t138;
t133 = t84 * qJD(3);
t115 = t85 * t133;
t116 = t81 * t140;
t39 = -qJD(2) * t115 + qJD(4) * t116 - t84 * t130;
t3 = t40 * pkin(4) + t39 * qJ(5) - t60 * qJD(5) + t23;
t174 = t3 * t81;
t173 = t3 * t84;
t135 = qJD(4) * t85;
t171 = -qJ(5) * t139 + t85 * qJD(5) - (-t82 * t133 - t81 * t135) * pkin(7) - t176;
t35 = -qJD(3) * pkin(3) - t44;
t58 = t116 - t133;
t11 = t58 * pkin(4) - t60 * qJ(5) + t35;
t170 = t11 * t60;
t169 = t23 * t81;
t168 = t23 * t84;
t167 = t39 * t81;
t166 = t58 * t72;
t165 = t60 * t58;
t164 = t60 * t72;
t163 = t72 * t81;
t162 = t72 * t84;
t161 = t79 * t83;
t160 = t79 * t86;
t88 = qJD(2) ^ 2;
t159 = t79 * t88;
t157 = t81 * t85;
t156 = t84 * t68;
t155 = t84 * t85;
t87 = qJD(3) ^ 2;
t153 = t87 * t82;
t152 = t87 * t85;
t151 = -pkin(4) * t139 + (-t82 * t134 + t84 * t135) * pkin(7) + t177;
t99 = pkin(4) * t81 - qJ(5) * t84;
t150 = t81 * qJD(5) + t72 * t99 + t45;
t62 = t103 * qJD(2);
t149 = t84 * t44 + t81 * t62;
t147 = pkin(7) * t155 + t81 * t68;
t77 = t82 ^ 2;
t146 = -t85 ^ 2 + t77;
t145 = qJD(2) * pkin(2);
t9 = -t81 * t36 + t84 * t46;
t144 = qJD(5) - t9;
t129 = pkin(8) * t163;
t128 = pkin(8) * t162;
t127 = t83 * t159;
t126 = pkin(8) * t139;
t125 = pkin(8) * t133;
t124 = t83 * t141;
t122 = t72 * t137;
t121 = t82 * t137;
t120 = t72 * t136;
t110 = t82 * t123;
t109 = t58 * t117;
t108 = t60 * t117;
t107 = t85 * t123;
t104 = qJ(5) * t112;
t66 = -t117 - t145;
t102 = -t66 - t117;
t6 = t72 * pkin(4) + t144;
t101 = t6 * t84 - t7 * t81;
t100 = t84 * pkin(4) + t81 * qJ(5);
t98 = -t81 * t44 + t84 * t62;
t97 = qJD(2) * t77 - t72 * t85;
t96 = pkin(7) + t99;
t51 = t85 * t161 + t158;
t28 = t84 * t160 + t51 * t81;
t29 = -t81 * t160 + t51 * t84;
t50 = t82 * t161 - t80 * t85;
t94 = -t10 * t72 + t111;
t93 = -t46 * t136 + t36 * t137 - t84 * t22 - t81 * t43;
t91 = qJD(3) * (-t102 - t145);
t27 = qJD(3) * t51 + t110;
t26 = -qJD(3) * t50 + t107;
t4 = qJD(4) * t29 - t84 * t124 + t26 * t81;
t90 = -t28 * t112 + t27 * t58 + t4 * t72 + t50 * t40;
t5 = -qJD(4) * t28 + t81 * t124 + t26 * t84;
t89 = t29 * t112 - t27 * t60 + t50 * t39 - t5 * t72;
t67 = -pkin(3) - t100;
t47 = t96 * t82;
t38 = -t156 + (pkin(7) * t81 + pkin(4)) * t85;
t37 = -t85 * qJ(5) + t147;
t25 = t60 * pkin(4) + t58 * qJ(5);
t16 = -t39 - t166;
t15 = (t100 * qJD(4) - qJD(5) * t84) * t82 + t96 * t138;
t14 = -pkin(4) * t140 - t98;
t13 = qJ(5) * t140 + t149;
t1 = -t72 * qJD(5) + t104 - t93;
t8 = [0, 0, -t127, -t86 * t159, 0, 0, 0, 0, 0, -t85 * t127 + (-t27 - t110) * qJD(3), t82 * t127 + (-t26 - t107) * qJD(3), 0, 0, 0, 0, 0, t90, -t89, t90, -t28 * t39 - t29 * t40 + t4 * t60 - t5 * t58, t89, t1 * t29 + t11 * t27 + t2 * t28 + t3 * t50 + t6 * t4 + t7 * t5; 0, 0, 0, 0, 0.2e1 * t85 * t112, -0.2e1 * t146 * t131, t152, -t153, 0, -pkin(7) * t152 + t82 * t91, pkin(7) * t153 + t85 * t91, -t39 * t84 * t82 + (t115 - t121) * t60, (-t58 * t84 - t60 * t81) * t138 + (t167 - t40 * t84 + (t58 * t81 - t60 * t84) * qJD(4)) * t82, t72 * t121 + t39 * t85 + (t60 * t82 + t84 * t97) * qJD(3), t72 * t119 + t40 * t85 + (-t58 * t82 - t81 * t97) * qJD(3), (-t72 - t132) * t139, t177 * t72 + (t35 * t134 + (qJD(3) * t58 + t120) * pkin(7) - t111) * t85 + (-t109 + t35 * t136 + pkin(7) * t40 + t169 + (-pkin(7) * t163 + (-pkin(7) * t157 + t156) * qJD(2) + t9) * qJD(3)) * t82, t176 * t72 + (t35 * t133 + (qJD(3) * t60 - t122) * pkin(7) - t93) * t85 + (-t108 - t35 * t137 - pkin(7) * t39 + t168 + (-pkin(7) * t162 - t147 * qJD(2) - t10) * qJD(3)) * t82, t15 * t58 + t47 * t40 + (t11 * t134 + t2) * t85 + t151 * t72 + (-t109 + t11 * t136 + t174 + (-qJD(2) * t38 - t6) * qJD(3)) * t82, -t37 * t40 - t38 * t39 + t151 * t60 + t171 * t58 + t101 * t138 + (-t1 * t81 + t2 * t84 + (-t6 * t81 - t7 * t84) * qJD(4)) * t82, -t15 * t60 + t47 * t39 + (-t11 * t133 - t1) * t85 + t171 * t72 + (t108 + t11 * t137 - t173 + (qJD(2) * t37 + t7) * qJD(3)) * t82, t1 * t37 + t2 * t38 + t3 * t47 - t171 * t7 + t151 * t6 + (t15 - t106) * t11; 0, 0, 0, 0, -t82 * t88 * t85, t146 * t88, 0, 0, 0, t45 * qJD(3) - t66 * t140 - t23, t102 * t132, -t60 * t162 - t167, (-t39 + t166) * t84 + (-t40 + t164) * t81, -t120 + (t72 * t155 + (-t60 + t134) * t82) * qJD(2), t122 + (-t72 * t157 + (t58 + t133) * t82) * qJD(2), t72 * t140, -pkin(3) * t40 - t168 + t98 * t72 - t45 * t58 + (t35 * t81 + t128) * qJD(4) + (-t9 * t82 + (-t35 * t85 - t126) * t81) * qJD(2), pkin(3) * t39 + t169 - t149 * t72 - t45 * t60 + (t35 * t84 - t129) * qJD(4) + (-t35 * t155 + (t10 - t125) * t82) * qJD(2), -t14 * t72 - t173 + t67 * t40 - t150 * t58 + (t11 * t81 + t128) * qJD(4) + (t6 * t82 + (-t11 * t85 - t126) * t81) * qJD(2), t13 * t58 - t14 * t60 + (t1 - t72 * t6 + (qJD(4) * t60 - t40) * pkin(8)) * t84 + ((qJD(4) * t58 - t39) * pkin(8) + t178) * t81, t13 * t72 - t174 + t67 * t39 + t150 * t60 + (-t11 * t84 + t129) * qJD(4) + (t11 * t155 + (-t7 + t125) * t82) * qJD(2), -t7 * t13 - t6 * t14 + t3 * t67 - t150 * t11 + (t101 * qJD(4) + t1 * t84 + t2 * t81) * pkin(8); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t165, -t58 ^ 2 + t175, t16, -t164 - t40, t112, -t35 * t60 + t94, t35 * t58 - t9 * t72 + t93, -t25 * t58 + 0.2e1 * t105 - t170 + t94, pkin(4) * t39 - t40 * qJ(5) + (-t10 + t7) * t60 + (t6 - t144) * t58, 0.2e1 * t104 - t11 * t58 + t25 * t60 + (-0.2e1 * qJD(5) + t9) * t72 - t93, -t2 * pkin(4) + t1 * qJ(5) - t6 * t10 - t11 * t25 + t144 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t112 + t165, t16, -t72 ^ 2 - t175, t170 + t178;];
tauc_reg = t8;
