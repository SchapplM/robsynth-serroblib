% Calculate minimal parameter regressor of coriolis joint torque vector for
% S6RPRPRP7
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
% tauc_reg [6x24]
%   minimal parameter regressor of coriolis joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 03:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc_reg = S6RPRPRP7_coriolisvecJ_fixb_regmin_slag_vp(qJ, qJD, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRP7_coriolisvecJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From coriolisvec_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:22:46
% EndTime: 2019-03-09 03:22:51
% DurationCPUTime: 1.62s
% Computational Cost: add. (3002->272), mult. (6467->381), div. (0->0), fcn. (4311->6), ass. (0->146)
t116 = cos(qJ(5));
t115 = sin(qJ(3));
t117 = cos(qJ(3));
t175 = sin(pkin(9));
t176 = cos(pkin(9));
t124 = t115 * t176 + t117 * t175;
t76 = t124 * qJD(1);
t197 = qJD(5) + t76;
t141 = t116 * t197;
t114 = sin(qJ(5));
t142 = qJD(3) * t175;
t129 = qJD(1) * t142;
t143 = qJD(3) * t176;
t130 = qJD(1) * t143;
t70 = t115 * t129 - t117 * t130;
t179 = t114 * t70;
t199 = t141 * t197 - t179;
t84 = -t115 * t175 + t117 * t176;
t163 = qJD(1) * t115;
t118 = -pkin(1) - pkin(7);
t93 = qJD(1) * t118 + qJD(2);
t73 = -qJ(4) * t163 + t115 * t93;
t147 = t176 * t73;
t162 = qJD(1) * t117;
t74 = -qJ(4) * t162 + t117 * t93;
t69 = qJD(3) * pkin(3) + t74;
t32 = t175 * t69 + t147;
t25 = qJD(3) * pkin(8) + t32;
t79 = t84 * qJD(1);
t87 = pkin(3) * t163 + qJD(1) * qJ(2) + qJD(4);
t33 = pkin(4) * t76 - pkin(8) * t79 + t87;
t13 = t114 * t33 + t116 * t25;
t155 = t117 * qJD(4);
t161 = qJD(3) * t115;
t121 = -t93 * t161 + (qJ(4) * t161 - t155) * qJD(1);
t157 = t115 * qJD(4);
t160 = qJD(3) * t117;
t54 = t93 * t160 + (-qJ(4) * t160 - t157) * qJD(1);
t18 = t121 * t175 + t176 * t54;
t110 = qJD(1) * qJD(2);
t153 = qJD(1) * qJD(3);
t146 = t117 * t153;
t177 = pkin(3) * t146 + t110;
t181 = t115 * t130 + t117 * t129;
t22 = -pkin(4) * t70 + pkin(8) * t181 + t177;
t20 = t116 * t22;
t122 = -qJD(5) * t13 - t114 * t18 + t20;
t156 = t116 * qJD(3);
t159 = qJD(5) * t114;
t26 = -qJD(5) * t156 + t116 * t181 + t159 * t79;
t58 = qJD(3) * t114 + t116 * t79;
t1 = -pkin(5) * t70 + qJ(6) * t26 - qJD(6) * t58 + t122;
t56 = t114 * t79 - t156;
t8 = -qJ(6) * t56 + t13;
t196 = t197 * t8 + t1;
t195 = t58 ^ 2;
t194 = 0.2e1 * t110;
t12 = -t114 * t25 + t116 * t33;
t7 = -qJ(6) * t58 + t12;
t5 = pkin(5) * t197 + t7;
t193 = t5 - t7;
t102 = pkin(3) * t175 + pkin(8);
t167 = qJ(6) + t102;
t139 = qJD(5) * t167;
t174 = qJ(6) * t114;
t66 = t175 * t73;
t41 = t176 * t74 - t66;
t43 = pkin(3) * t162 + pkin(4) * t79 + pkin(8) * t76;
t183 = t114 * t43 + t116 * t41;
t192 = t116 * qJD(6) - t114 * t139 - t174 * t76 - t183;
t173 = qJ(6) * t116;
t39 = t116 * t43;
t191 = -pkin(5) * t79 - t116 * t139 - t76 * t173 - t39 + (-qJD(6) + t41) * t114;
t17 = -t121 * t176 + t175 * t54;
t190 = t17 * t84;
t168 = qJ(4) - t118;
t88 = t168 * t115;
t89 = t168 * t117;
t52 = -t175 * t89 - t176 * t88;
t189 = t52 * t70;
t188 = t56 * t76;
t187 = t56 * t79;
t186 = t58 * t79;
t185 = t124 * t70;
t158 = qJD(5) * t116;
t148 = t114 * t181;
t27 = qJD(5) * t58 - t148;
t184 = -t114 * t27 - t158 * t56;
t169 = pkin(3) * t115 + qJ(2);
t48 = pkin(4) * t124 - pkin(8) * t84 + t169;
t49 = t116 * t52;
t182 = t114 * t48 + t49;
t180 = t114 * t58;
t64 = t116 * t70;
t178 = t26 * t114;
t119 = qJD(3) ^ 2;
t172 = t119 * t115;
t171 = t119 * t117;
t120 = qJD(1) ^ 2;
t170 = t120 * qJ(2);
t165 = t115 ^ 2 - t117 ^ 2;
t164 = -t119 - t120;
t154 = pkin(3) * t160 + qJD(2);
t71 = t161 * t168 - t155;
t72 = -qJD(3) * t89 - t157;
t37 = t175 * t71 + t176 * t72;
t77 = t115 * t142 - t117 * t143;
t78 = -t115 * t143 - t117 * t142;
t42 = -pkin(4) * t77 - pkin(8) * t78 + t154;
t152 = t114 * t42 + t116 * t37 + t158 * t48;
t151 = 0.2e1 * qJD(1);
t150 = t84 * t159;
t149 = t84 * t158;
t140 = t197 * t114;
t138 = -qJD(5) * t124 - qJD(1);
t103 = -pkin(3) * t176 - pkin(4);
t137 = qJD(5) * t102 * t197 + t17;
t126 = t114 * t22 + t116 * t18 + t158 * t33 - t159 * t25;
t2 = -qJ(6) * t27 - qJD(6) * t56 + t126;
t136 = -t197 * t5 + t2;
t31 = t176 * t69 - t66;
t24 = -qJD(3) * pkin(4) - t31;
t135 = t24 * t78 + t190;
t134 = -t124 * t26 - t58 * t77;
t133 = -t26 * t84 + t58 * t78;
t132 = -t124 * t27 + t56 * t77;
t131 = -t197 * t78 + t70 * t84;
t36 = t175 * t72 - t176 * t71;
t40 = t175 * t74 + t147;
t51 = -t175 * t88 + t176 * t89;
t128 = -qJ(6) * t78 - qJD(6) * t84;
t127 = -t64 + (-t114 * t76 - t159) * t197;
t125 = t102 * t70 + t197 * t24;
t10 = pkin(5) * t27 + t17;
t123 = -t124 * t18 - t31 * t78 + t32 * t77 + t190;
t82 = t167 * t116;
t81 = t167 * t114;
t55 = t56 ^ 2;
t46 = t116 * t48;
t35 = t116 * t42;
t15 = pkin(5) * t56 + qJD(6) + t24;
t14 = -t174 * t84 + t182;
t11 = pkin(5) * t124 - t114 * t52 - t173 * t84 + t46;
t4 = -qJ(6) * t149 + (-qJD(5) * t52 + t128) * t114 + t152;
t3 = -pkin(5) * t77 - t114 * t37 + t35 + t128 * t116 + (-t49 + (qJ(6) * t84 - t48) * t114) * qJD(5);
t6 = [0, 0, 0, 0, t194, qJ(2) * t194, -0.2e1 * t115 * t146, 0.2e1 * t165 * t153, -t172, -t171, 0, -t118 * t172 + (qJ(2) * t160 + qJD(2) * t115) * t151, -t118 * t171 + (-qJ(2) * t161 + qJD(2) * t117) * t151, -t181 * t51 + t36 * t79 - t37 * t76 + t123 + t189, t154 * t87 + t169 * t177 + t17 * t51 + t18 * t52 - t31 * t36 + t32 * t37, t116 * t133 - t150 * t58 (-t116 * t56 - t180) * t78 + (t178 - t116 * t27 + (t114 * t56 - t116 * t58) * qJD(5)) * t84, -t116 * t131 - t150 * t197 + t134, t114 * t131 - t149 * t197 + t132, -t197 * t77 - t185 (-t52 * t158 + t35) * t197 - t46 * t70 + (-t25 * t158 + t20) * t124 - t12 * t77 + t36 * t56 + t51 * t27 + t24 * t149 + ((-qJD(5) * t48 - t37) * t197 + t189 + (-qJD(5) * t33 - t18) * t124 + t135) * t114 -(-t159 * t52 + t152) * t197 + t182 * t70 - t126 * t124 + t13 * t77 + t36 * t58 - t51 * t26 - t24 * t150 + t135 * t116, t11 * t26 - t14 * t27 - t3 * t58 - t4 * t56 + (-t114 * t8 - t116 * t5) * t78 + (-t1 * t116 - t114 * t2 + (t114 * t5 - t116 * t8) * qJD(5)) * t84, t2 * t14 + t8 * t4 + t1 * t11 + t5 * t3 + t10 * (pkin(5) * t114 * t84 + t51) + t15 * ((t114 * t78 + t149) * pkin(5) + t36); 0, 0, 0, 0, -t120, -t170, 0, 0, 0, 0, 0, t164 * t115, t164 * t117, t181 * t84 + t76 * t77 - t78 * t79 + t185, -qJD(1) * t87 - t123, 0, 0, 0, 0, 0, t124 * t179 - t27 * t84 - t56 * t78 + (t114 * t77 + t116 * t138) * t197, t124 * t64 + (-t114 * t138 + t116 * t77) * t197 - t133 (-t138 * t58 + t132) * t116 + (-t138 * t56 + t134) * t114, -t10 * t84 - t15 * t78 + (t124 * t2 + t138 * t5 - t77 * t8) * t116 + (-t1 * t124 + t138 * t8 + t5 * t77) * t114; 0, 0, 0, 0, 0, 0, t117 * t120 * t115, -t165 * t120, 0, 0, 0, -t117 * t170, t115 * t170 (t32 - t40) * t79 - (-t41 + t31) * t76 + (t175 * t70 + t176 * t181) * pkin(3), t31 * t40 - t32 * t41 + (-t162 * t87 - t17 * t176 + t175 * t18) * pkin(3), t141 * t58 - t178 (-t26 - t188) * t116 - t197 * t180 + t184, -t186 + t199, t127 + t187, -t197 * t79, t103 * t27 - t12 * t79 - t39 * t197 - t40 * t56 - t137 * t116 + (t197 * t41 + t125) * t114, -t103 * t26 + t114 * t137 + t116 * t125 + t13 * t79 + t183 * t197 - t40 * t58, -t114 * t196 + t116 * t136 - t191 * t58 - t192 * t56 - t26 * t81 - t27 * t82, t2 * t82 - t1 * t81 + t10 * (-pkin(5) * t116 + t103) + t192 * t8 + t191 * t5 + (pkin(5) * t140 - t40) * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76 ^ 2 - t79 ^ 2, t31 * t79 + t32 * t76 + t177, 0, 0, 0, 0, 0, t127 - t187, -t186 - t199 (t26 - t188) * t116 + t58 * t140 + t184, t114 * t136 + t116 * t196 - t15 * t79; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t58 * t56, -t55 + t195, t197 * t56 - t26, t148 + (-qJD(5) + t197) * t58, -t70, t13 * t197 - t24 * t58 + t122, t12 * t197 + t24 * t56 - t126, pkin(5) * t26 - t193 * t56, t193 * t8 + (-t15 * t58 + t1) * pkin(5); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t55 - t195, t5 * t58 + t56 * t8 + t10;];
tauc_reg  = t6;
