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
% Datum: 2022-01-20 09:49
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:49:04
% EndTime: 2022-01-20 09:49:08
% DurationCPUTime: 1.31s
% Computational Cost: add. (5980->208), mult. (8272->297), div. (0->0), fcn. (5095->10), ass. (0->141)
t140 = sin(qJ(5));
t131 = qJDD(4) + qJDD(5);
t134 = qJD(1) + qJD(3);
t144 = cos(qJ(5));
t145 = cos(qJ(4));
t141 = sin(qJ(4));
t177 = t134 * t141;
t95 = -t144 * t145 * t134 + t140 * t177;
t97 = (t145 * t140 + t141 * t144) * t134;
t75 = t97 * t95;
t189 = -t75 + t131;
t192 = t140 * t189;
t191 = t144 * t189;
t137 = -g(3) + qJDD(2);
t130 = t134 ^ 2;
t132 = qJDD(1) + qJDD(3);
t142 = sin(qJ(3));
t146 = cos(qJ(3));
t138 = sin(pkin(9));
t139 = cos(pkin(9));
t149 = qJD(1) ^ 2;
t143 = sin(qJ(1));
t147 = cos(qJ(1));
t160 = t147 * g(1) + t143 * g(2);
t154 = -t149 * pkin(1) - t160;
t165 = t143 * g(1) - t147 * g(2);
t180 = qJDD(1) * pkin(1);
t155 = t165 + t180;
t151 = -t138 * t154 + t139 * t155;
t150 = qJDD(1) * pkin(2) + t151;
t173 = t138 * t155 + t139 * t154;
t79 = -t149 * pkin(2) + t173;
t63 = t142 * t150 + t146 * t79;
t51 = -t130 * pkin(3) + t132 * pkin(7) + t63;
t43 = -t145 * t137 + t141 * t51;
t44 = t141 * t137 + t145 * t51;
t24 = t141 * t43 + t145 * t44;
t172 = qJD(4) * t134;
t167 = t145 * t172;
t175 = t141 * t132;
t104 = t167 + t175;
t122 = t145 * t132;
t168 = t141 * t172;
t159 = t122 - t168;
t66 = -t95 * qJD(5) + t144 * t104 + t140 * t159;
t133 = qJD(4) + qJD(5);
t92 = t133 * t95;
t190 = -t92 + t66;
t188 = -t43 + (-t104 + t167) * pkin(8);
t93 = t95 ^ 2;
t94 = t97 ^ 2;
t129 = t133 ^ 2;
t119 = t145 * t130 * t141;
t171 = qJDD(4) + t119;
t152 = t171 * pkin(4) + t188;
t115 = qJD(4) * pkin(4) - pkin(8) * t177;
t136 = t145 ^ 2;
t124 = t136 * t130;
t37 = -pkin(4) * t124 + t159 * pkin(8) - qJD(4) * t115 + t44;
t18 = t140 * t37 - t144 * t152;
t183 = t144 * t37;
t19 = t140 * t152 + t183;
t7 = t140 * t19 - t144 * t18;
t187 = t141 * t7;
t62 = -t142 * t79 + t146 * t150;
t50 = -t132 * pkin(3) - t130 * pkin(7) - t62;
t186 = -pkin(3) * t50 + pkin(7) * t24;
t39 = -t159 * pkin(4) - pkin(8) * t124 + t115 * t177 + t50;
t185 = t140 * t39;
t72 = t75 + t131;
t184 = t140 * t72;
t182 = t144 * t39;
t181 = t144 * t72;
t179 = t133 * t140;
t178 = t133 * t144;
t176 = t141 * t171;
t114 = qJDD(4) - t119;
t174 = t145 * t114;
t103 = 0.2e1 * t167 + t175;
t135 = t141 ^ 2;
t123 = t135 * t130;
t148 = qJD(4) ^ 2;
t116 = -t123 - t148;
t87 = -t141 * t116 - t174;
t170 = -pkin(3) * t103 + pkin(7) * t87 + t141 * t50;
t105 = t122 - 0.2e1 * t168;
t117 = -t124 - t148;
t86 = t145 * t117 - t176;
t169 = pkin(3) * t105 + pkin(7) * t86 - t145 * t50;
t162 = t140 * t104 - t144 * t159;
t153 = (-qJD(5) + t133) * t97 - t162;
t59 = t92 + t66;
t29 = t140 * t153 - t144 * t59;
t30 = t140 * t59 + t144 * t153;
t12 = -t141 * t29 + t145 * t30;
t67 = -t93 - t94;
t8 = t140 * t18 + t144 * t19;
t166 = t141 * (-pkin(8) * t29 - t7) + t145 * (-pkin(4) * t67 + pkin(8) * t30 + t8) - pkin(3) * t67 + pkin(7) * t12;
t70 = -t129 - t93;
t45 = t140 * t70 + t191;
t46 = t144 * t70 - t192;
t26 = -t141 * t45 + t145 * t46;
t54 = (qJD(5) + t133) * t97 + t162;
t164 = t141 * (-pkin(8) * t45 + t185) + t145 * (-pkin(4) * t54 + pkin(8) * t46 - t182) - pkin(3) * t54 + pkin(7) * t26;
t88 = -t94 - t129;
t60 = t144 * t88 - t184;
t61 = -t140 * t88 - t181;
t32 = -t141 * t60 + t145 * t61;
t163 = t141 * (-pkin(8) * t60 + t182) + t145 * (-pkin(4) * t190 + pkin(8) * t61 + t185) - pkin(3) * t190 + pkin(7) * t32;
t109 = (t135 + t136) * t132;
t112 = t123 + t124;
t161 = pkin(3) * t112 + pkin(7) * t109 + t24;
t110 = -t146 * t130 - t142 * t132;
t157 = t142 * t130 - t146 * t132;
t3 = t145 * t8 - t187;
t156 = pkin(7) * t3 - pkin(8) * t187 - pkin(3) * t39 + t145 * (-pkin(4) * t39 + pkin(8) * t8);
t90 = -t94 + t129;
t89 = t93 - t129;
t85 = t176 + t145 * (-t123 + t148);
t84 = t141 * (t124 - t148) + t174;
t81 = (t104 + t167) * t141;
t80 = t105 * t145;
t78 = t142 * t109 + t146 * t112;
t76 = t145 * t103 + t141 * t105;
t74 = t94 - t93;
t69 = -t146 * t103 + t142 * t87;
t68 = t146 * t105 + t142 * t86;
t65 = -t97 * qJD(5) - t162;
t40 = (t141 * (t140 * t97 - t144 * t95) + t145 * (-t140 * t95 - t144 * t97)) * t133;
t35 = t141 * (t144 * t89 - t184) + t145 * (t140 * t89 + t181);
t34 = t141 * (-t140 * t90 + t191) + t145 * (t144 * t90 + t192);
t33 = t142 * t63 + t146 * t62;
t28 = t141 * (t144 * t66 - t97 * t179) + t145 * (t140 * t66 + t97 * t178);
t27 = t141 * (-t140 * t65 + t95 * t178) + t145 * (t144 * t65 + t95 * t179);
t20 = t142 * t32 - t146 * t190;
t16 = t142 * t26 - t146 * t54;
t15 = t142 * t24 - t146 * t50;
t11 = t141 * (-t140 * t190 - t144 * t54) + t145 * (-t140 * t54 + t144 * t190);
t9 = t142 * t12 - t146 * t67;
t1 = t142 * t3 - t146 * t39;
t2 = [0, 0, 0, 0, 0, qJDD(1), t165, t160, 0, 0, 0, 0, 0, 0, 0, qJDD(1), t138 * t160 + (t165 + 0.2e1 * t180) * t139, pkin(1) * (-t138 * qJDD(1) - t139 * t149) - t173, 0, pkin(1) * (t138 * t173 + t139 * t151), 0, 0, 0, 0, 0, t132, pkin(1) * (t138 * t110 - t139 * t157) - pkin(2) * t157 + t62, pkin(1) * (t139 * t110 + t138 * t157) + pkin(2) * t110 - t63, 0, pkin(1) * (t138 * (-t142 * t62 + t146 * t63) + t139 * t33) + pkin(2) * t33, t81, t76, t85, t80, t84, 0, pkin(1) * (t138 * (-t142 * t105 + t146 * t86) + t139 * t68) + pkin(2) * t68 + t169, pkin(1) * (t138 * (t142 * t103 + t146 * t87) + t139 * t69) + pkin(2) * t69 + t170, pkin(1) * (t138 * (t146 * t109 - t142 * t112) + t139 * t78) + pkin(2) * t78 + t161, pkin(1) * (t138 * (t142 * t50 + t146 * t24) + t139 * t15) + pkin(2) * t15 + t186, t28, t11, t34, t27, t35, t40, pkin(1) * (t138 * (t142 * t54 + t146 * t26) + t139 * t16) + pkin(2) * t16 + t164, pkin(1) * (t138 * (t142 * t190 + t146 * t32) + t139 * t20) + pkin(2) * t20 + t163, pkin(1) * (t138 * (t146 * t12 + t142 * t67) + t139 * t9) + pkin(2) * t9 + t166, pkin(1) * (t138 * (t142 * t39 + t146 * t3) + t139 * t1) + pkin(2) * t1 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, 0, 0, 0, 0, 0, 0, 0, 0, t137, 0, 0, 0, 0, 0, 0, t141 * t117 + t145 * t171, -t141 * t114 + t145 * t116, 0, t141 * t44 - t145 * t43, 0, 0, 0, 0, 0, 0, t141 * t46 + t145 * t45, t141 * t61 + t145 * t60, t141 * t30 + t145 * t29, t141 * t8 + t145 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t62, -t63, 0, 0, t81, t76, t85, t80, t84, 0, t169, t170, t161, t186, t28, t11, t34, t27, t35, t40, t164, t163, t166, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t123 - t124, t175, t119, t122, qJDD(4), -t43, -t44, 0, 0, t75, t74, t59, -t75, t153, t131, pkin(4) * t45 - t18, -t183 - t140 * t188 + (-t140 * t171 + t60) * pkin(4), pkin(4) * t29, pkin(4) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t59, -t75, t153, t131, -t18, -t19, 0, 0;];
tauJ_reg = t2;
