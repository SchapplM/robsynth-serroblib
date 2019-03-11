% Calculate minimal parameter regressor of inverse dynamics joint torque vector for
% S6RPPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3]';
% 
% Output:
% tau_reg [6x26]
%   minimal parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 01:36
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau_reg = S6RPPPRR4_invdynJ_fixb_regmin_slag_vp(qJ, qJD, qJDD, g, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_regmin_slag_vp: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_regmin_slag_vp: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR4_invdynJ_fixb_regmin_slag_vp: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR4_invdynJ_fixb_regmin_slag_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPPRR4_invdynJ_fixb_regmin_slag_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_joint_fixb_regressor_minpar_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:35:55
% EndTime: 2019-03-09 01:35:58
% DurationCPUTime: 1.35s
% Computational Cost: add. (1232->268), mult. (2086->330), div. (0->0), fcn. (1327->8), ass. (0->154)
t83 = cos(qJ(5));
t148 = qJD(6) * t83;
t191 = qJD(1) * t148 + qJDD(5);
t133 = qJD(1) * qJD(5);
t81 = sin(qJ(5));
t125 = t81 * t133;
t136 = t83 * qJDD(1);
t190 = -qJD(5) * qJD(6) - t125 + t136;
t189 = -qJDD(1) * qJ(4) - qJD(1) * qJD(4);
t188 = qJDD(1) * pkin(3) + qJDD(4);
t141 = qJ(2) * qJD(1);
t84 = -pkin(1) - pkin(2);
t49 = t84 * qJD(1) + qJD(2);
t77 = sin(pkin(9));
t78 = cos(pkin(9));
t25 = -t77 * t141 + t78 * t49;
t22 = qJD(1) * pkin(3) + qJD(4) - t25;
t19 = qJD(1) * pkin(7) + t22;
t13 = -t81 * qJD(3) + t83 * t19;
t147 = t13 * qJD(5);
t135 = qJ(2) * qJDD(1);
t48 = t84 * qJDD(1) + qJDD(2);
t162 = t77 * t135 - t78 * t48;
t134 = qJD(1) * qJD(2);
t53 = t77 * t134;
t20 = -t53 - t162;
t18 = -t20 + t188;
t15 = qJDD(1) * pkin(7) + t18;
t171 = t77 * t49;
t72 = qJD(1) * qJ(4);
t112 = -t81 * pkin(5) + t83 * pkin(8);
t157 = t78 * qJ(2);
t99 = t112 + t157;
t16 = t99 * qJD(1) + t171 - t72;
t120 = -qJDD(5) * pkin(8) - qJD(6) * t16 - t83 * qJDD(3) - t81 * t15 - t147;
t150 = qJD(5) * t83;
t154 = qJD(2) * t77;
t176 = sin(qJ(1));
t177 = cos(qJ(1));
t31 = -t176 * t78 + t177 * t77;
t181 = g(1) * t31;
t14 = t83 * qJD(3) + t81 * t19;
t146 = t14 * qJD(5);
t2 = -qJDD(5) * pkin(5) + t81 * qJDD(3) - t83 * t15 + t146;
t170 = t77 * t84;
t24 = -qJ(4) + t99 + t170;
t39 = -t77 * qJ(2) + t78 * t84;
t35 = pkin(3) - t39;
t28 = pkin(7) + t35;
t122 = t83 * t133;
t137 = t81 * qJDD(1);
t29 = -qJDD(6) + t122 + t137;
t51 = t81 * qJD(1) - qJD(6);
t7 = -qJD(5) * pkin(5) - t13;
t187 = (qJD(6) * t24 + t28 * t150) * t51 + (t7 * qJD(5) + t51 * t154 + t28 * t29 - t120) * t81 - t181 - t2 * t83;
t138 = qJDD(5) * t83;
t85 = qJD(5) ^ 2;
t86 = qJD(1) ^ 2;
t158 = t85 + t86;
t186 = t158 * t81 - t138;
t139 = qJDD(5) * t81;
t185 = t158 * t83 + t139;
t26 = t78 * t141 + t171;
t23 = t26 - t72;
t40 = t157 + t170;
t32 = -qJ(4) + t40;
t184 = (qJD(1) * t32 - t154 + t23) * qJD(5) - qJDD(5) * t28;
t30 = -t176 * t77 - t177 * t78;
t111 = -g(2) * t30 + t181;
t113 = -pkin(5) * t83 - pkin(8) * t81;
t183 = t111 * t83 - (pkin(8) * qJD(6) + t113 * qJD(1)) * t51 + g(3) * t81 + t2;
t182 = g(1) * t30;
t180 = g(2) * t31;
t80 = sin(qJ(6));
t82 = cos(qJ(6));
t9 = -t190 * t82 + t191 * t80;
t179 = t83 * t9;
t178 = t9 * t80;
t175 = t30 * t81;
t142 = t82 * qJD(5);
t155 = qJD(1) * t83;
t33 = t80 * t155 + t142;
t174 = t33 * t51;
t143 = t80 * qJD(5);
t34 = t82 * t155 - t143;
t173 = t34 * t51;
t172 = t34 * t82;
t169 = t80 * t29;
t168 = t80 * t81;
t167 = t81 * t82;
t166 = t82 * t29;
t165 = t83 * t33;
t164 = t83 * t34;
t163 = t78 * t135 + t77 * t48;
t161 = t177 * pkin(1) + t176 * qJ(2);
t160 = g(1) * t176 - g(2) * t177;
t74 = t83 ^ 2;
t159 = t81 ^ 2 - t74;
t156 = pkin(1) * qJDD(1);
t153 = qJD(5) * t33;
t152 = qJD(5) * t34;
t151 = qJD(5) * t81;
t149 = qJD(6) * t51;
t145 = t23 * qJD(1);
t144 = t78 * qJD(2);
t75 = qJDD(3) + g(3);
t131 = t177 * pkin(2) + t161;
t130 = t51 * t143;
t129 = t51 * t142;
t128 = 0.2e1 * t134;
t8 = qJD(5) * pkin(8) + t14;
t126 = t28 * t51 + t8;
t123 = t78 * t134;
t119 = -t163 - t189;
t118 = -0.2e1 * t122;
t117 = qJD(6) * t81 - qJD(1);
t116 = qJDD(2) - t156;
t115 = t9 - t129;
t114 = -t176 * pkin(1) + t177 * qJ(2);
t110 = t180 + t182;
t10 = t190 * t80 + t191 * t82;
t109 = -t10 - t130;
t108 = -qJD(6) * t8 + t180;
t91 = t113 * qJD(5) + t144;
t107 = -t24 * t29 - (-qJD(4) + t91) * t51;
t106 = t25 * t77 - t26 * t78;
t105 = -qJD(5) * t19 - t75;
t21 = t123 + t163;
t104 = t82 * t149 + t169;
t103 = t80 * t149 - t166;
t100 = -t111 + t162;
t98 = g(1) * t177 + g(2) * t176;
t97 = -t32 * qJDD(1) - t110;
t96 = -t176 * pkin(2) + t114;
t95 = t100 + t188;
t94 = t104 - t153;
t93 = t103 + t152;
t92 = -g(2) * t175 - g(3) * t83 + t120;
t90 = pkin(8) * t29 + (-t13 - t7) * t51;
t89 = qJD(3) * qJD(5) + t111 - t145 - t15;
t17 = t21 + t189;
t50 = -qJD(4) + t144;
t88 = -t50 * qJD(1) - t28 * t85 - t17 + t97;
t44 = t85 * t81 - t138;
t43 = t85 * t83 + t139;
t42 = t78 * qJDD(1) + t77 * t86;
t41 = t77 * qJDD(1) - t78 * t86;
t12 = t31 * t167 - t30 * t80;
t11 = -t31 * t168 - t30 * t82;
t6 = t91 * qJD(1) + t112 * qJDD(1) - t119;
t5 = t82 * t6;
t4 = t80 * t16 + t82 * t8;
t3 = t82 * t16 - t80 * t8;
t1 = [qJDD(1), t160, t98, -qJDD(2) + 0.2e1 * t156 + t160, t128 - t98 + 0.2e1 * t135, -t116 * pkin(1) - g(1) * t114 - g(2) * t161 + (t128 + t135) * qJ(2), -t39 * qJDD(1) + t100 + 0.2e1 * t53, t40 * qJDD(1) + t110 + 0.2e1 * t123 + t163, -g(1) * t96 - g(2) * t131 - t106 * qJD(2) + t20 * t39 + t21 * t40, -t35 * qJDD(1) - 0.2e1 * t53 - t95 (-t50 - t144) * qJD(1) + t97 + t119, t17 * t32 + t23 * t50 + t18 * t35 + t22 * t154 - g(1) * (t31 * pkin(3) + t30 * qJ(4) + t96) - g(2) * (-t30 * pkin(3) + t31 * qJ(4) + t131) t74 * qJDD(1) + t81 * t118, 0.2e1 * t159 * t133 - 0.2e1 * t81 * t136, t44, t43, 0, -t184 * t83 + t88 * t81, t184 * t81 + t88 * t83, -t82 * t179 + (-t81 * t142 - t80 * t148) * t34 (t33 * t82 + t34 * t80) * t151 + (-t10 * t82 + t178 + (t33 * t80 - t172) * qJD(6)) * t83 (-t9 - t129) * t81 + (-t103 + t152) * t83 (-t10 + t130) * t81 + (-t104 - t153) * t83, t51 * t150 + t29 * t81, -g(2) * t12 + (-t28 * t153 - t5) * t81 + (-t3 * qJD(5) + t28 * t10 + t33 * t154) * t83 + (-g(1) * t175 + (t126 * t81 - t7 * t83) * qJD(6) + t107) * t82 + t187 * t80, -t28 * t34 * t151 - g(2) * t11 + (t4 * qJD(5) + t34 * t154 - t28 * t9) * t83 + (t7 * t148 + (-t126 * qJD(6) + t182 + t6) * t81 - t107) * t80 + t187 * t82; 0, 0, 0, -qJDD(1), -t86, -t86 * qJ(2) + t116 - t160, -t42, t41, t106 * qJD(1) + t20 * t78 + t21 * t77 - t160, t42, -t41, t17 * t77 - t18 * t78 + (-t22 * t77 - t23 * t78) * qJD(1) - t160, 0, 0, 0, 0, 0 (t118 - t137) * t77 + t186 * t78 (0.2e1 * t125 - t136) * t77 + t185 * t78, 0, 0, 0, 0, 0, t103 * t77 + (t109 * t83 - t94 * t81) * t78 + ((-t77 * t168 + t78 * t82) * t51 - t77 * t165) * qJD(1), t104 * t77 + (t115 * t83 + t93 * t81) * t78 + (-(t77 * t167 + t78 * t80) * t51 - t77 * t164) * qJD(1); 0, 0, 0, 0, 0, 0, 0, 0, t75, 0, 0, t75, 0, 0, 0, 0, 0, -t43, t44, 0, 0, 0, 0, 0, t109 * t81 + t94 * t83, t115 * t81 - t83 * t93; 0, 0, 0, 0, 0, 0, 0, 0, 0, -qJDD(1), -t86, t53 + t95 + t145, 0, 0, 0, 0, 0, -t186, -t185, 0, 0, 0, 0, 0, t83 * t10 + (-t153 + t169) * t81 + (t117 * t82 + t83 * t143) * t51, -t179 + (-t152 + t166) * t81 + (-t117 * t80 + t83 * t142) * t51; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t83 * t86 * t81, -t159 * t86, -t136, t137, qJDD(5), t105 * t81 - t89 * t83 + t146, t105 * t83 + t89 * t81 + t147, t51 * t172 + t178 (t9 - t174) * t82 + (t10 - t173) * t80 (t51 * t167 - t164) * qJD(1) - t104 (-t51 * t168 + t165) * qJD(1) + t103, -t51 * t155, pkin(5) * t10 + t14 * t33 + t3 * t155 - t183 * t82 + t90 * t80, -pkin(5) * t9 + t14 * t34 - t4 * t155 + t183 * t80 + t90 * t82; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t34 * t33, -t33 ^ 2 + t34 ^ 2, t9 + t174, t10 + t173, -t29, -g(1) * t11 + t108 * t82 + t7 * t34 - t4 * t51 + t92 * t80 + t5, g(1) * t12 - t3 * t51 - t7 * t33 + (-t108 - t6) * t80 + t92 * t82;];
tau_reg  = t1;
