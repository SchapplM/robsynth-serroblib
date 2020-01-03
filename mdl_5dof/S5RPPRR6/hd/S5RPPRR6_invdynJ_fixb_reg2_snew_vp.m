% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRR6_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR6_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:58:08
% EndTime: 2019-12-31 17:58:14
% DurationCPUTime: 2.14s
% Computational Cost: add. (6012->276), mult. (13459->403), div. (0->0), fcn. (9347->10), ass. (0->162)
t135 = sin(pkin(9));
t130 = t135 ^ 2;
t137 = cos(pkin(9));
t131 = t137 ^ 2;
t193 = qJD(1) ^ 2;
t118 = (t130 + t131) * t193;
t140 = sin(qJ(5));
t141 = sin(qJ(4));
t143 = cos(qJ(4));
t155 = t135 * t143 + t137 * t141;
t112 = t155 * qJD(1);
t142 = cos(qJ(5));
t97 = -t142 * qJD(4) + t140 * t112;
t99 = t140 * qJD(4) + t142 * t112;
t78 = t99 * t97;
t167 = t137 * qJDD(1);
t168 = t135 * qJDD(1);
t157 = t141 * t168 - t143 * t167;
t172 = t112 * qJD(4);
t91 = -t157 - t172;
t84 = qJDD(5) - t91;
t198 = -t78 + t84;
t204 = t140 * t198;
t174 = qJD(1) * t137;
t177 = t135 * t141;
t110 = qJD(1) * t177 - t143 * t174;
t94 = t112 * t110;
t195 = qJDD(4) - t94;
t203 = t141 * t195;
t202 = t142 * t198;
t201 = t143 * t195;
t136 = sin(pkin(8));
t190 = sin(qJ(1));
t191 = cos(qJ(1));
t152 = t190 * g(1) - t191 * g(2);
t153 = t191 * g(1) + t190 * g(2);
t116 = -t193 * pkin(1) - t153;
t138 = cos(pkin(8));
t176 = t138 * t116;
t150 = -t136 * t152 - t176;
t175 = -g(3) + qJDD(2);
t160 = t137 * t175;
t162 = t136 * pkin(1) + qJ(3);
t199 = t162 + pkin(6);
t147 = t160 + (-t199 * qJDD(1) + (-(2 * qJD(3)) + (t137 * pkin(3) + pkin(2)) * qJD(1)) * qJD(1) + t150) * t135;
t171 = t131 * t193;
t151 = qJDD(1) * pkin(1) + t152;
t149 = qJDD(1) * qJ(3) + t136 * t151 + t176;
t192 = 2 * qJD(3);
t74 = t137 * (-t193 * pkin(2) + t149) + t135 * t175 + t174 * t192;
t67 = -pkin(3) * t171 + pkin(6) * t167 + t74;
t41 = t141 * t147 + t143 * t67;
t164 = pkin(1) * t138 + pkin(2);
t113 = t138 * t151;
t159 = -t136 * t116 + t113;
t83 = -qJDD(1) * pkin(2) - t193 * qJ(3) + qJDD(3) - t159;
t200 = -t164 * qJDD(1) + t162 * t118 + t83;
t106 = qJD(5) + t110;
t109 = t155 * qJDD(1);
t173 = t110 * qJD(4);
t93 = t109 - t173;
t161 = -t142 * qJDD(4) + t140 * t93;
t48 = (qJD(5) - t106) * t99 + t161;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t105 = t106 ^ 2;
t107 = t110 ^ 2;
t108 = t112 ^ 2;
t40 = t141 * t67 - t143 * t147;
t20 = t141 * t41 - t143 * t40;
t189 = t135 * t20;
t144 = qJD(4) ^ 2;
t85 = t110 * pkin(4) - t112 * pkin(7);
t29 = -qJDD(4) * pkin(4) - t144 * pkin(7) + t112 * t85 + t40;
t188 = t140 * t29;
t56 = t78 + t84;
t187 = t140 * t56;
t72 = -pkin(3) * t167 + t83 + (-t130 * t193 - t171) * pkin(6);
t186 = t141 * t72;
t88 = qJDD(4) + t94;
t185 = t141 * t88;
t184 = t142 * t29;
t183 = t142 * t56;
t181 = t143 * t72;
t180 = t143 * t88;
t179 = t106 * t140;
t178 = t106 * t142;
t169 = qJD(5) + t106;
t166 = t141 * t78;
t165 = t143 * t78;
t163 = -pkin(4) * t143 - pkin(3);
t30 = -t144 * pkin(4) + qJDD(4) * pkin(7) - t110 * t85 + t41;
t32 = (-t93 + t173) * pkin(7) + (-t91 + t172) * pkin(4) + t72;
t14 = t140 * t30 - t142 * t32;
t15 = t140 * t32 + t142 * t30;
t6 = t140 * t14 + t142 * t15;
t73 = -t160 + ((-pkin(2) * qJD(1) + t192) * qJD(1) + t149) * t135;
t45 = t135 * t73 + t137 * t74;
t21 = t141 * t40 + t143 * t41;
t5 = -t142 * t14 + t140 * t15;
t154 = -t140 * qJDD(4) - t142 * t93;
t65 = -t97 * qJD(5) - t154;
t127 = t131 * qJDD(1);
t126 = t130 * qJDD(1);
t117 = t127 + t126;
t102 = -t108 - t144;
t101 = -t108 + t144;
t100 = t107 - t144;
t92 = t109 - 0.2e1 * t173;
t90 = t157 + 0.2e1 * t172;
t86 = -t144 - t107;
t81 = t106 * t97;
t80 = -t96 + t105;
t79 = t95 - t105;
t77 = -t107 - t108;
t76 = t96 - t95;
t70 = -t96 - t105;
t69 = -t141 * t102 - t180;
t68 = t143 * t102 - t185;
t66 = -t105 - t95;
t64 = -t99 * qJD(5) - t161;
t62 = t95 + t96;
t61 = t141 * t109 - t143 * t157;
t60 = -t143 * t109 - t141 * t157;
t59 = t143 * t86 - t203;
t58 = t141 * t86 + t201;
t54 = (t140 * t99 - t142 * t97) * t106;
t53 = t169 * t97 + t154;
t52 = t65 + t81;
t51 = t65 - t81;
t49 = -t169 * t99 - t161;
t47 = t142 * t65 - t99 * t179;
t46 = -t140 * t64 + t97 * t178;
t44 = -t135 * t68 + t137 * t69;
t43 = t142 * t79 - t187;
t42 = -t140 * t80 + t202;
t38 = -t140 * t70 - t183;
t37 = t142 * t70 - t187;
t36 = -t135 * t60 + t137 * t61;
t35 = t142 * t66 - t204;
t34 = t140 * t66 + t202;
t33 = -t135 * t58 + t137 * t59;
t28 = t140 * t52 - t142 * t48;
t27 = -t140 * t51 + t142 * t49;
t26 = -t140 * t48 - t142 * t52;
t25 = -t141 * t53 + t143 * t38;
t24 = t141 * t38 + t143 * t53;
t23 = -t141 * t49 + t143 * t35;
t22 = t141 * t35 + t143 * t49;
t19 = -t141 * t62 + t143 * t28;
t18 = t141 * t28 + t143 * t62;
t17 = -pkin(7) * t37 + t184;
t16 = -pkin(7) * t34 + t188;
t12 = -pkin(4) * t37 + t15;
t11 = -pkin(4) * t34 + t14;
t10 = -t135 * t24 + t137 * t25;
t9 = -t135 * t22 + t137 * t23;
t8 = t137 * t21 - t189;
t4 = t141 * t29 + t143 * t6;
t3 = t141 * t6 - t143 * t29;
t2 = -pkin(7) * t26 - t5;
t1 = [0, 0, 0, 0, 0, qJDD(1), t152, t153, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (t138 * qJDD(1) - t136 * t193) + t159, (-0.2e1 * t136 * qJDD(1) - t138 * t193) * pkin(1) + t150, 0, pkin(1) * (t136 ^ 2 * t151 + t138 * t113), t126, 0.2e1 * t135 * t167, 0, t127, 0, 0, -t200 * t137, t200 * t135, pkin(2) * t118 + qJ(3) * t117 + pkin(1) * (t136 * t117 + t138 * t118) + t45, -pkin(2) * t83 + qJ(3) * t45 + pkin(1) * (t136 * t45 - t138 * t83), t135 * (-t141 * t172 + t143 * t93) + t137 * (t141 * t93 + t143 * t172), t135 * (-t141 * t92 - t143 * t90) + t137 * (-t141 * t90 + t143 * t92), t135 * (-t141 * t101 + t201) + t137 * (t143 * t101 + t203), t135 * (-t141 * t91 + t143 * t173) + t137 * (t141 * t173 + t143 * t91), t135 * (t143 * t100 - t185) + t137 * (t141 * t100 + t180), (t135 * (-t110 * t143 + t112 * t141) + t137 * (-t110 * t141 - t112 * t143)) * qJD(4), t135 * (-pkin(6) * t58 + t186) + t137 * (-pkin(3) * t90 + pkin(6) * t59 - t181) - pkin(2) * t90 + qJ(3) * t33 + pkin(1) * (t136 * t33 - t138 * t90), t135 * (-pkin(6) * t68 + t181) + t137 * (-pkin(3) * t92 + pkin(6) * t69 + t186) - pkin(2) * t92 + qJ(3) * t44 + pkin(1) * (t136 * t44 - t138 * t92), t135 * (-pkin(6) * t60 - t20) + t137 * (-pkin(3) * t77 + pkin(6) * t61 + t21) - pkin(2) * t77 + qJ(3) * t36 + pkin(1) * (t136 * t36 - t138 * t77), -pkin(6) * t189 + t137 * (-pkin(3) * t72 + pkin(6) * t21) - pkin(2) * t72 + qJ(3) * t8 + pkin(1) * (t136 * t8 - t138 * t72), t135 * (t143 * t47 + t166) + t137 * (t141 * t47 - t165), t135 * (t141 * t76 + t143 * t27) + t137 * (t141 * t27 - t143 * t76), t135 * (t141 * t52 + t143 * t42) + t137 * (t141 * t42 - t143 * t52), t135 * (t143 * t46 - t166) + t137 * (t141 * t46 + t165), t135 * (-t141 * t48 + t143 * t43) + t137 * (t141 * t43 + t143 * t48), t135 * (t141 * t84 + t143 * t54) + t137 * (t141 * t54 - t143 * t84), t135 * (-pkin(6) * t22 - t141 * t11 + t143 * t16) + t137 * (-pkin(3) * t34 + pkin(6) * t23 + t143 * t11 + t141 * t16) - pkin(2) * t34 + qJ(3) * t9 + pkin(1) * (t136 * t9 - t138 * t34), t135 * (-pkin(6) * t24 - t141 * t12 + t143 * t17) + t137 * (-pkin(3) * t37 + pkin(6) * t25 + t143 * t12 + t141 * t17) - pkin(2) * t37 + qJ(3) * t10 + pkin(1) * (t136 * t10 - t138 * t37), t135 * (-pkin(6) * t18 + t143 * t2) + t137 * (pkin(6) * t19 + t141 * t2) + t162 * (-t135 * t18 + t137 * t19) + (pkin(4) * t177 + t137 * t163 - t164) * t26, (t135 * (pkin(4) * t141 - pkin(7) * t143) + t137 * (-pkin(7) * t141 + t163) - t164) * t5 + t199 * (-t135 * t3 + t137 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t175, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135 * t74 - t137 * t73, 0, 0, 0, 0, 0, 0, t135 * t59 + t137 * t58, t135 * t69 + t137 * t68, t135 * t61 + t137 * t60, t135 * t21 + t137 * t20, 0, 0, 0, 0, 0, 0, t135 * t23 + t137 * t22, t135 * t25 + t137 * t24, t135 * t19 + t137 * t18, t135 * t4 + t137 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t167, t168, -t118, t83, 0, 0, 0, 0, 0, 0, t90, t92, t77, t72, 0, 0, 0, 0, 0, 0, t34, t37, t26, t5; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t94, t108 - t107, t109, -t94, -t157, qJDD(4), -t40, -t41, 0, 0, t140 * t65 + t99 * t178, t140 * t49 + t142 * t51, t142 * t80 + t204, t142 * t64 + t97 * t179, t140 * t79 + t183, (-t140 * t97 - t142 * t99) * t106, pkin(4) * t49 + pkin(7) * t35 - t184, pkin(4) * t53 + pkin(7) * t38 + t188, pkin(4) * t62 + pkin(7) * t28 + t6, -pkin(4) * t29 + pkin(7) * t6; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t78, t76, t52, -t78, -t48, t84, -t14, -t15, 0, 0;];
tauJ_reg = t1;
