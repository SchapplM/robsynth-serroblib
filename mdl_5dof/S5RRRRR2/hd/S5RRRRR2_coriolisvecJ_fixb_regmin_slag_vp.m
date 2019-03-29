% Calculate minimal parameter regressor of coriolis joint torque vector for
% S5RRRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
% 
% Output:
% tauc_reg [5x27]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-29 15:26
% Revision: 932832b1be1be80f59b7f1a581a1a8f328bdb39d (2019-03-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S5RRRRR2_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [2x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-29 15:26:31
% EndTime: 2019-03-29 15:26:35
% DurationCPUTime: 1.54s
% Computational Cost: add. (1847->181), mult. (3878->299), div. (0->0), fcn. (2870->8), ass. (0->132)
t175 = cos(qJ(4));
t98 = cos(qJ(3));
t134 = t175 * t98;
t94 = sin(qJ(4));
t95 = sin(qJ(3));
t107 = -t94 * t95 + t134;
t90 = qJD(1) + qJD(2);
t179 = t107 * t90;
t89 = qJD(3) + qJD(4);
t35 = t179 * t89;
t69 = t175 * t95 + t94 * t98;
t178 = t89 * t69;
t36 = t178 * t90;
t132 = qJD(4) * t175;
t151 = pkin(1) * qJD(1);
t139 = t98 * t151;
t99 = cos(qJ(2));
t148 = qJD(2) * t99;
t96 = sin(qJ(2));
t140 = t96 * t151;
t128 = t95 * t140;
t115 = qJD(3) * t128;
t157 = t96 * t98;
t58 = (-qJD(3) * t157 - t95 * t148) * t151;
t153 = -t94 * t115 - t175 * t58;
t73 = qJD(3) * pkin(2) - t128;
t160 = t94 * t73;
t22 = (t96 * t132 + t94 * t148) * t139 + qJD(4) * t160 + t153;
t64 = t69 * t90;
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t116 = -t97 * t64 - t93 * t89;
t17 = -t116 * qJD(5) + t93 * t35;
t114 = t134 * t151;
t76 = t96 * t114;
t50 = t76 + t160;
t163 = t90 * t98;
t67 = -pkin(2) * t163 - t99 * t151;
t119 = t93 * t50 - t97 * t67;
t129 = t94 * t139;
t122 = t96 * t129;
t21 = -qJD(4) * t122 + t114 * t148 - t175 * t115 + t73 * t132 + t94 * t58;
t150 = pkin(1) * qJD(2);
t138 = t96 * t150;
t127 = qJD(1) * t138;
t147 = qJD(3) * t95;
t137 = t90 * t147;
t65 = pkin(2) * t137 + t127;
t177 = -t119 * qJD(5) + t97 * t21 + t93 * t65;
t176 = pkin(1) * t96;
t144 = qJD(5) * t97;
t145 = qJD(5) * t93;
t16 = t89 * t144 - t64 * t145 + t97 * t35;
t174 = t16 * t93;
t42 = t93 * t64 - t97 * t89;
t56 = qJD(5) - t179;
t173 = t42 * t56;
t172 = t116 * t56;
t45 = t89 * t107;
t171 = t45 * t93;
t170 = t45 * t97;
t169 = t56 * t64;
t168 = t64 * t179;
t167 = t67 * t64;
t166 = t69 * t93;
t165 = t69 * t97;
t164 = t90 * t95;
t161 = t93 * t36;
t156 = t97 * t36;
t155 = t67 * t45 + t65 * t69;
t154 = -t107 * t65 + t178 * t67;
t152 = t95 ^ 2 - t98 ^ 2;
t100 = qJD(3) ^ 2;
t149 = t100 * t95;
t87 = t100 * t98;
t146 = qJD(3) * t99;
t143 = qJD(5) * t98;
t142 = -qJD(1) - t90;
t136 = t95 * t146;
t49 = -t175 * t73 + t122;
t40 = t49 * t145;
t41 = t49 * t144;
t133 = t49 * (-t56 - t179);
t131 = t56 * t97;
t31 = t97 * t50 + t93 * t67;
t130 = t107 * t177 + t22 * t165 + t49 * t170 - t178 * t31;
t126 = t22 * t93 + t31 * t64 + t41;
t103 = -t31 * qJD(5) - t93 * t21 + t97 * t65;
t125 = -t103 * t107 - t119 * t178 + t22 * t166 + t49 * t171 + t69 * t41;
t124 = (-qJD(2) + t90) * t151;
t123 = t142 * t150;
t54 = t69 * t140;
t120 = -t179 * t49 - t54 * t56;
t61 = t107 * t176;
t83 = -t99 * pkin(1) - t98 * pkin(2);
t118 = t97 * t61 + t93 * t83;
t117 = -t93 * t61 + t97 * t83;
t113 = qJD(5) * t94 + t164;
t111 = t96 * t124;
t110 = t99 * t124;
t109 = t119 * t64 - t22 * t97 + t40;
t108 = -t69 * t145 + t170;
t106 = t99 * t69;
t105 = t99 * t107;
t101 = -t179 * t67 - t21;
t88 = t90 ^ 2;
t79 = t95 * t127;
t75 = 0.2e1 * t98 * t137;
t74 = pkin(2) * t147 + t138;
t60 = t69 * t176;
t59 = -0.2e1 * t152 * t90 * qJD(3);
t55 = t105 * t151;
t53 = t106 * t151;
t52 = t94 * t128 - t76;
t39 = t178 * t89;
t38 = t45 * t89;
t27 = -t179 ^ 2 + t64 ^ 2;
t26 = (qJD(2) * t106 + t45 * t96) * pkin(1);
t25 = (qJD(2) * t105 - t178 * t96) * pkin(1);
t24 = t64 * t89 - t36;
t13 = t35 * t69 + t64 * t45;
t12 = -t107 * t36 + t178 * t56;
t9 = t116 * t64 + t56 * t131 + t161;
t8 = -t56 ^ 2 * t93 + t42 * t64 + t156;
t7 = -t116 * t131 + t174;
t6 = -t108 * t116 + t16 * t165;
t5 = t107 * t35 - t178 * t64 + t179 * t45 - t69 * t36;
t4 = -t69 * t161 + t17 * t107 - t42 * t178 + (-t69 * t144 - t171) * t56;
t3 = -t107 * t16 + t108 * t56 - t116 * t178 + t69 * t156;
t2 = (t16 - t173) * t97 + (-t17 + t172) * t93;
t1 = (t116 * t93 - t42 * t97) * t45 + (-t174 - t17 * t97 + (t116 * t97 + t42 * t93) * qJD(5)) * t69;
t10 = [0, 0, 0, 0, t96 * t123, t99 * t123, t75, t59, t87, -t149, 0 (-t96 * t87 + t142 * t136 + (t142 * t157 - t136) * qJD(2)) * pkin(1), t79 + ((qJD(2) * t90 + t100) * t96 * t95 - 0.2e1 * t146 * t163) * pkin(1), t13, t5, t38, -t39, 0, -t179 * t74 - t26 * t89 + t83 * t36 + t154, -t25 * t89 + t83 * t35 + t74 * t64 + t155, t6, t1, t3, t4, t12 (-qJD(5) * t118 - t93 * t25 + t97 * t74) * t56 + t117 * t36 + t26 * t42 + t60 * t17 + t125 -(t97 * t25 + t93 * t74) * t56 - t118 * t36 - t26 * t116 + t60 * t16 + (-t117 * t56 - t49 * t166) * qJD(5) + t130; 0, 0, 0, 0, t111, t110, t75, t59, t87, -t149, 0, t98 * t111, -t90 * t128 + t79, t13, t5, t38, -t39, 0, t179 * t140 + t53 * t89 + (-t147 * t179 - t36 * t98) * pkin(2) + t154, -t64 * t140 + t55 * t89 + (t64 * t147 - t35 * t98) * pkin(2) + t155, t6, t1, t3, t4, t12 -(t140 * t97 - t93 * t55) * t56 - t53 * t42 + (-t98 * t156 + (t143 * t93 + t147 * t97) * t56) * pkin(2) + t125, -t69 * t40 + (t140 * t93 + t97 * t55) * t56 + t53 * t116 + (t98 * t161 + (t143 * t97 - t147 * t93) * t56) * pkin(2) + t130; 0, 0, 0, 0, 0, 0, -t95 * t88 * t98, t152 * t88, 0, 0, 0, t95 * t110, t98 * t110, -t168, t27, 0, t24, 0, -t129 * t148 + pkin(2) * t179 * t164 - t52 * t89 - t167 + (-t76 + (-pkin(2) * t89 - t73) * t94) * qJD(4) - t153, -t54 * t89 + (-t89 * t132 - t64 * t164) * pkin(2) + t101, t7, t2, t9, t8, -t169, t52 * t42 + t120 * t93 + (-t175 * t17 + (qJD(4) * t42 - t161) * t94 + (-t113 * t97 - t132 * t93) * t56) * pkin(2) + t109, -t52 * t116 + t120 * t97 + (-t175 * t16 + (-qJD(4) * t116 - t156) * t94 + (t113 * t93 - t132 * t97) * t56) * pkin(2) + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t168, t27, 0, t24, 0, t50 * t89 - t167 - t22, -t49 * t89 + t101, t7, t2, t9, t8, -t169, t133 * t93 - t50 * t42 + t109, t116 * t50 + t133 * t97 + t126; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t116 * t42, t116 ^ 2 - t42 ^ 2, t16 + t173, -t17 - t172, t36, t116 * t49 + t31 * t56 + t103, -t119 * t56 + t49 * t42 - t177;];
tauc_reg  = t10;
