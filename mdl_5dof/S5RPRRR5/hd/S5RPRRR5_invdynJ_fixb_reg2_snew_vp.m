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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:16:30
% EndTime: 2019-12-05 18:16:36
% DurationCPUTime: 1.23s
% Computational Cost: add. (5980->209), mult. (8272->297), div. (0->0), fcn. (5095->10), ass. (0->141)
t141 = sin(qJ(5));
t132 = qJDD(4) + qJDD(5);
t135 = qJD(1) + qJD(3);
t145 = cos(qJ(5));
t146 = cos(qJ(4));
t142 = sin(qJ(4));
t178 = t135 * t142;
t95 = -t145 * t146 * t135 + t141 * t178;
t97 = (t146 * t141 + t142 * t145) * t135;
t75 = t97 * t95;
t190 = -t75 + t132;
t193 = t141 * t190;
t192 = t145 * t190;
t138 = -g(1) + qJDD(2);
t131 = t135 ^ 2;
t133 = qJDD(1) + qJDD(3);
t143 = sin(qJ(3));
t147 = cos(qJ(3));
t139 = sin(pkin(9));
t140 = cos(pkin(9));
t150 = qJD(1) ^ 2;
t144 = sin(qJ(1));
t148 = cos(qJ(1));
t161 = t144 * g(2) - t148 * g(3);
t155 = -t150 * pkin(1) + t161;
t173 = t148 * g(2) + t144 * g(3);
t181 = qJDD(1) * pkin(1);
t160 = t173 + t181;
t153 = -t139 * t155 + t140 * t160;
t151 = qJDD(1) * pkin(2) + t153;
t174 = t139 * t160 + t140 * t155;
t79 = -t150 * pkin(2) + t174;
t63 = t143 * t151 + t147 * t79;
t51 = -t131 * pkin(3) + t133 * pkin(7) + t63;
t43 = -t146 * t138 + t142 * t51;
t44 = t142 * t138 + t146 * t51;
t24 = t142 * t43 + t146 * t44;
t172 = qJD(4) * t135;
t167 = t146 * t172;
t176 = t142 * t133;
t104 = t167 + t176;
t122 = t146 * t133;
t168 = t142 * t172;
t159 = t122 - t168;
t66 = -t95 * qJD(5) + t145 * t104 + t141 * t159;
t134 = qJD(4) + qJD(5);
t92 = t134 * t95;
t191 = -t92 + t66;
t189 = -t43 + (-t104 + t167) * pkin(8);
t93 = t95 ^ 2;
t94 = t97 ^ 2;
t130 = t134 ^ 2;
t119 = t146 * t131 * t142;
t171 = qJDD(4) + t119;
t152 = t171 * pkin(4) + t189;
t115 = qJD(4) * pkin(4) - pkin(8) * t178;
t137 = t146 ^ 2;
t124 = t137 * t131;
t37 = -pkin(4) * t124 + t159 * pkin(8) - qJD(4) * t115 + t44;
t18 = t141 * t37 - t145 * t152;
t184 = t145 * t37;
t19 = t141 * t152 + t184;
t7 = t141 * t19 - t145 * t18;
t188 = t142 * t7;
t62 = -t143 * t79 + t147 * t151;
t50 = -t133 * pkin(3) - t131 * pkin(7) - t62;
t187 = -pkin(3) * t50 + pkin(7) * t24;
t39 = -t159 * pkin(4) - pkin(8) * t124 + t115 * t178 + t50;
t186 = t141 * t39;
t72 = t75 + t132;
t185 = t141 * t72;
t183 = t145 * t39;
t182 = t145 * t72;
t180 = t134 * t141;
t179 = t134 * t145;
t177 = t142 * t171;
t114 = qJDD(4) - t119;
t175 = t146 * t114;
t103 = 0.2e1 * t167 + t176;
t136 = t142 ^ 2;
t123 = t136 * t131;
t149 = qJD(4) ^ 2;
t116 = -t123 - t149;
t87 = -t142 * t116 - t175;
t170 = -pkin(3) * t103 + pkin(7) * t87 + t142 * t50;
t105 = t122 - 0.2e1 * t168;
t117 = -t124 - t149;
t86 = t146 * t117 - t177;
t169 = pkin(3) * t105 + pkin(7) * t86 - t146 * t50;
t163 = t141 * t104 - t145 * t159;
t154 = (-qJD(5) + t134) * t97 - t163;
t59 = t92 + t66;
t29 = t141 * t154 - t145 * t59;
t30 = t141 * t59 + t145 * t154;
t12 = -t142 * t29 + t146 * t30;
t67 = -t93 - t94;
t8 = t141 * t18 + t145 * t19;
t166 = t142 * (-pkin(8) * t29 - t7) + t146 * (-pkin(4) * t67 + pkin(8) * t30 + t8) - pkin(3) * t67 + pkin(7) * t12;
t70 = -t130 - t93;
t45 = t141 * t70 + t192;
t46 = t145 * t70 - t193;
t26 = -t142 * t45 + t146 * t46;
t54 = (qJD(5) + t134) * t97 + t163;
t165 = t142 * (-pkin(8) * t45 + t186) + t146 * (-pkin(4) * t54 + pkin(8) * t46 - t183) - pkin(3) * t54 + pkin(7) * t26;
t88 = -t94 - t130;
t60 = t145 * t88 - t185;
t61 = -t141 * t88 - t182;
t32 = -t142 * t60 + t146 * t61;
t164 = t142 * (-pkin(8) * t60 + t183) + t146 * (-pkin(4) * t191 + pkin(8) * t61 + t186) - pkin(3) * t191 + pkin(7) * t32;
t109 = (t136 + t137) * t133;
t112 = t123 + t124;
t162 = pkin(3) * t112 + pkin(7) * t109 + t24;
t110 = -t147 * t131 - t143 * t133;
t157 = t143 * t131 - t147 * t133;
t3 = t146 * t8 - t188;
t156 = pkin(7) * t3 - pkin(8) * t188 - pkin(3) * t39 + t146 * (-pkin(4) * t39 + pkin(8) * t8);
t90 = -t94 + t130;
t89 = t93 - t130;
t85 = t177 + t146 * (-t123 + t149);
t84 = t142 * (t124 - t149) + t175;
t81 = (t104 + t167) * t142;
t80 = t105 * t146;
t78 = t143 * t109 + t147 * t112;
t76 = t146 * t103 + t142 * t105;
t74 = t94 - t93;
t69 = -t147 * t103 + t143 * t87;
t68 = t147 * t105 + t143 * t86;
t65 = -t97 * qJD(5) - t163;
t40 = (t142 * (t141 * t97 - t145 * t95) + t146 * (-t141 * t95 - t145 * t97)) * t134;
t35 = t142 * (t145 * t89 - t185) + t146 * (t141 * t89 + t182);
t34 = t142 * (-t141 * t90 + t192) + t146 * (t145 * t90 + t193);
t33 = t143 * t63 + t147 * t62;
t28 = t142 * (t145 * t66 - t97 * t180) + t146 * (t141 * t66 + t97 * t179);
t27 = t142 * (-t141 * t65 + t95 * t179) + t146 * (t145 * t65 + t95 * t180);
t20 = t143 * t32 - t147 * t191;
t16 = t143 * t26 - t147 * t54;
t15 = t143 * t24 - t147 * t50;
t11 = t142 * (-t141 * t191 - t145 * t54) + t146 * (-t141 * t54 + t145 * t191);
t9 = t143 * t12 - t147 * t67;
t1 = t143 * t3 - t147 * t39;
t2 = [0, 0, 0, 0, 0, qJDD(1), t173, -t161, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t139 * t161 + (t173 + 0.2e1 * t181) * t140, pkin(1) * (-t139 * qJDD(1) - t140 * t150) - t174, 0, pkin(1) * (t139 * t174 + t140 * t153), 0, 0, 0, 0, 0, t133, pkin(1) * (t139 * t110 - t140 * t157) - pkin(2) * t157 + t62, pkin(1) * (t140 * t110 + t139 * t157) + pkin(2) * t110 - t63, 0, pkin(1) * (t139 * (-t143 * t62 + t147 * t63) + t140 * t33) + pkin(2) * t33, t81, t76, t85, t80, t84, 0, pkin(1) * (t139 * (-t143 * t105 + t147 * t86) + t140 * t68) + pkin(2) * t68 + t169, pkin(1) * (t139 * (t143 * t103 + t147 * t87) + t140 * t69) + pkin(2) * t69 + t170, pkin(1) * (t139 * (t147 * t109 - t143 * t112) + t140 * t78) + pkin(2) * t78 + t162, pkin(1) * (t139 * (t143 * t50 + t147 * t24) + t140 * t15) + pkin(2) * t15 + t187, t28, t11, t34, t27, t35, t40, pkin(1) * (t139 * (t143 * t54 + t147 * t26) + t140 * t16) + pkin(2) * t16 + t165, pkin(1) * (t139 * (t143 * t191 + t147 * t32) + t140 * t20) + pkin(2) * t20 + t164, pkin(1) * (t139 * (t147 * t12 + t143 * t67) + t140 * t9) + pkin(2) * t9 + t166, pkin(1) * (t139 * (t143 * t39 + t147 * t3) + t140 * t1) + pkin(2) * t1 + t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, 0, 0, 0, t138, 0, 0, 0, 0, 0, 0, t142 * t117 + t146 * t171, -t142 * t114 + t146 * t116, 0, t142 * t44 - t146 * t43, 0, 0, 0, 0, 0, 0, t142 * t46 + t146 * t45, t142 * t61 + t146 * t60, t142 * t30 + t146 * t29, t142 * t8 + t146 * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, t62, -t63, 0, 0, t81, t76, t85, t80, t84, 0, t169, t170, t162, t187, t28, t11, t34, t27, t35, t40, t165, t164, t166, t156; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t119, t123 - t124, t176, t119, t122, qJDD(4), -t43, -t44, 0, 0, t75, t74, t59, -t75, t154, t132, pkin(4) * t45 - t18, -t184 - t141 * t189 + (-t141 * t171 + t60) * pkin(4), pkin(4) * t29, pkin(4) * t7; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t75, t74, t59, -t75, t154, t132, -t18, -t19, 0, 0;];
tauJ_reg = t2;
