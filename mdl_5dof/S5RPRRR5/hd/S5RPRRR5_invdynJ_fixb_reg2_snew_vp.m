% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRR5
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRR5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:54:09
% EndTime: 2020-01-03 11:54:14
% DurationCPUTime: 1.27s
% Computational Cost: add. (5980->208), mult. (8272->297), div. (0->0), fcn. (5095->10), ass. (0->141)
t139 = sin(qJ(5));
t130 = qJDD(4) + qJDD(5);
t133 = qJD(1) + qJD(3);
t143 = cos(qJ(5));
t144 = cos(qJ(4));
t140 = sin(qJ(4));
t176 = t133 * t140;
t95 = -t143 * t144 * t133 + t139 * t176;
t97 = (t144 * t139 + t140 * t143) * t133;
t75 = t97 * t95;
t188 = -t75 + t130;
t191 = t139 * t188;
t190 = t143 * t188;
t136 = -g(1) + qJDD(2);
t129 = t133 ^ 2;
t131 = qJDD(1) + qJDD(3);
t141 = sin(qJ(3));
t145 = cos(qJ(3));
t137 = sin(pkin(9));
t138 = cos(pkin(9));
t142 = sin(qJ(1));
t146 = cos(qJ(1));
t160 = -t146 * g(2) - t142 * g(3);
t179 = qJDD(1) * pkin(1);
t153 = t160 + t179;
t148 = qJD(1) ^ 2;
t159 = t142 * g(2) - t146 * g(3);
t154 = -t148 * pkin(1) - t159;
t150 = -t137 * t154 + t138 * t153;
t149 = qJDD(1) * pkin(2) + t150;
t172 = t137 * t153 + t138 * t154;
t79 = -t148 * pkin(2) + t172;
t63 = t141 * t149 + t145 * t79;
t51 = -t129 * pkin(3) + t131 * pkin(7) + t63;
t43 = -t144 * t136 + t140 * t51;
t44 = t140 * t136 + t144 * t51;
t24 = t140 * t43 + t144 * t44;
t171 = qJD(4) * t133;
t166 = t144 * t171;
t174 = t140 * t131;
t104 = t166 + t174;
t122 = t144 * t131;
t167 = t140 * t171;
t158 = t122 - t167;
t66 = -t95 * qJD(5) + t143 * t104 + t139 * t158;
t132 = qJD(4) + qJD(5);
t92 = t132 * t95;
t189 = -t92 + t66;
t187 = -t43 + (-t104 + t166) * pkin(8);
t93 = t95 ^ 2;
t94 = t97 ^ 2;
t128 = t132 ^ 2;
t119 = t144 * t129 * t140;
t170 = qJDD(4) + t119;
t151 = t170 * pkin(4) + t187;
t115 = qJD(4) * pkin(4) - pkin(8) * t176;
t135 = t144 ^ 2;
t124 = t135 * t129;
t37 = -pkin(4) * t124 + t158 * pkin(8) - qJD(4) * t115 + t44;
t18 = t139 * t37 - t143 * t151;
t182 = t143 * t37;
t19 = t139 * t151 + t182;
t7 = t139 * t19 - t143 * t18;
t186 = t140 * t7;
t62 = -t141 * t79 + t145 * t149;
t50 = -t131 * pkin(3) - t129 * pkin(7) - t62;
t185 = -pkin(3) * t50 + pkin(7) * t24;
t39 = -t158 * pkin(4) - pkin(8) * t124 + t115 * t176 + t50;
t184 = t139 * t39;
t72 = t75 + t130;
t183 = t139 * t72;
t181 = t143 * t39;
t180 = t143 * t72;
t178 = t132 * t139;
t177 = t132 * t143;
t175 = t140 * t170;
t114 = qJDD(4) - t119;
t173 = t144 * t114;
t103 = 0.2e1 * t166 + t174;
t134 = t140 ^ 2;
t123 = t134 * t129;
t147 = qJD(4) ^ 2;
t116 = -t123 - t147;
t87 = -t140 * t116 - t173;
t169 = -pkin(3) * t103 + pkin(7) * t87 + t140 * t50;
t105 = t122 - 0.2e1 * t167;
t117 = -t124 - t147;
t85 = t144 * t117 - t175;
t168 = pkin(3) * t105 + pkin(7) * t85 - t144 * t50;
t162 = t139 * t104 - t143 * t158;
t152 = (-qJD(5) + t132) * t97 - t162;
t59 = t92 + t66;
t29 = t139 * t152 - t143 * t59;
t30 = t139 * t59 + t143 * t152;
t12 = -t140 * t29 + t144 * t30;
t67 = -t93 - t94;
t8 = t139 * t18 + t143 * t19;
t165 = t140 * (-pkin(8) * t29 - t7) + t144 * (-pkin(4) * t67 + pkin(8) * t30 + t8) - pkin(3) * t67 + pkin(7) * t12;
t70 = -t128 - t93;
t45 = t139 * t70 + t190;
t46 = t143 * t70 - t191;
t26 = -t140 * t45 + t144 * t46;
t54 = (qJD(5) + t132) * t97 + t162;
t164 = t140 * (-pkin(8) * t45 + t184) + t144 * (-pkin(4) * t54 + pkin(8) * t46 - t181) - pkin(3) * t54 + pkin(7) * t26;
t88 = -t94 - t128;
t60 = t143 * t88 - t183;
t61 = -t139 * t88 - t180;
t32 = -t140 * t60 + t144 * t61;
t163 = t140 * (-pkin(8) * t60 + t181) + t144 * (-pkin(4) * t189 + pkin(8) * t61 + t184) - pkin(3) * t189 + pkin(7) * t32;
t109 = (t134 + t135) * t131;
t112 = t123 + t124;
t161 = pkin(3) * t112 + pkin(7) * t109 + t24;
t110 = -t145 * t129 - t141 * t131;
t156 = t141 * t129 - t145 * t131;
t3 = t144 * t8 - t186;
t155 = pkin(7) * t3 - pkin(8) * t186 - pkin(3) * t39 + t144 * (-pkin(4) * t39 + pkin(8) * t8);
t90 = -t94 + t128;
t89 = t93 - t128;
t86 = t173 + t140 * (t124 - t147);
t84 = t144 * (-t123 + t147) + t175;
t81 = t105 * t144;
t80 = (t104 + t166) * t140;
t78 = t141 * t109 + t145 * t112;
t76 = t144 * t103 + t140 * t105;
t74 = t94 - t93;
t69 = -t145 * t103 + t141 * t87;
t68 = t145 * t105 + t141 * t85;
t65 = -t97 * qJD(5) - t162;
t40 = (t144 * (-t139 * t95 - t143 * t97) + t140 * (t139 * t97 - t143 * t95)) * t132;
t35 = t144 * (t139 * t89 + t180) + t140 * (t143 * t89 - t183);
t34 = t144 * (t143 * t90 + t191) + t140 * (-t139 * t90 + t190);
t33 = t141 * t63 + t145 * t62;
t28 = t144 * (t139 * t66 + t97 * t177) + t140 * (t143 * t66 - t97 * t178);
t27 = t144 * (t143 * t65 + t95 * t178) + t140 * (-t139 * t65 + t95 * t177);
t20 = t141 * t32 - t145 * t189;
t16 = t141 * t26 - t145 * t54;
t15 = t141 * t24 - t145 * t50;
t11 = t144 * (-t139 * t54 + t143 * t189) + t140 * (-t139 * t189 - t143 * t54);
t9 = t141 * t12 - t145 * t67;
t1 = t141 * t3 - t145 * t39;
t2 = [0, 0, 0, 0, 0, qJDD(1), t160, t159, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t137 * t159 + (t160 + 0.2e1 * t179) * t138, pkin(1) * (-t137 * qJDD(1) - t138 * t148) - t172, 0, pkin(1) * (t137 * t172 + t138 * t150), 0, 0, 0, 0, 0, t131, pkin(1) * (t137 * t110 - t138 * t156) - pkin(2) * t156 + t62, pkin(1) * (t138 * t110 + t137 * t156) + pkin(2) * t110 - t63, 0, pkin(1) * (t137 * (-t141 * t62 + t145 * t63) + t138 * t33) + pkin(2) * t33, t80, t76, t84, t81, t86, 0, pkin(1) * (t137 * (-t141 * t105 + t145 * t85) + t138 * t68) + pkin(2) * t68 + t168, pkin(1) * (t137 * (t141 * t103 + t145 * t87) + t138 * t69) + pkin(2) * t69 + t169, pkin(1) * (t137 * (t145 * t109 - t141 * t112) + t138 * t78) + pkin(2) * t78 + t161, pkin(1) * (t137 * (t141 * t50 + t145 * t24) + t138 * t15) + pkin(2) * t15 + t185, t28, t11, t34, t27, t35, t40, pkin(1) * (t137 * (t141 * t54 + t145 * t26) + t138 * t16) + pkin(2) * t16 + t164, pkin(1) * (t137 * (t141 * t189 + t145 * t32) + t138 * t20) + pkin(2) * t20 + t163, pkin(1) * (t137 * (t145 * t12 + t141 * t67) + t138 * t9) + pkin(2) * t9 + t165, pkin(1) * (t137 * (t141 * t39 + t145 * t3) + t138 * t1) + pkin(2) * t1 + t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, 0, 0, 0, t136, 0, 0, 0, 0, 0, 0, t140 * t117 + t144 * t170, -t140 * t114 + t144 * t116, 0, t140 * t44 - t144 * t43, 0, 0, 0, 0, 0, 0, t140 * t46 + t144 * t45, t140 * t61 + t144 * t60, t140 * t30 + t144 * t29, t140 * t8 + t144 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t62, -t63, 0, 0, t80, t76, t84, t81, t86, 0, t168, t169, t161, t185, t28, t11, t34, t27, t35, t40, t164, t163, t165, t155; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t123 - t124, t174, t119, t122, qJDD(4), -t43, -t44, 0, 0, t75, t74, t59, -t75, t152, t130, pkin(4) * t45 - t18, -t182 - t139 * t187 + (-t139 * t170 + t60) * pkin(4), pkin(4) * t29, pkin(4) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t59, -t75, t152, t130, -t18, -t19, 0, 0;];
tauJ_reg = t2;
