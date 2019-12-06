% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5PRRPP2
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 16:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5PRRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:10:09
% EndTime: 2019-12-05 16:10:17
% DurationCPUTime: 1.90s
% Computational Cost: add. (2716->215), mult. (6114->267), div. (0->0), fcn. (4043->8), ass. (0->136)
t131 = sin(qJ(3));
t133 = cos(qJ(3));
t128 = cos(pkin(8));
t126 = sin(pkin(8));
t155 = qJD(2) * t131;
t96 = -t128 * t133 * qJD(2) + t126 * t155;
t160 = t131 * t128;
t98 = (t133 * t126 + t160) * qJD(2);
t175 = t98 * t96;
t181 = qJDD(3) + t175;
t170 = t126 * t181;
t135 = qJD(3) ^ 2;
t95 = t98 ^ 2;
t185 = -t95 - t135;
t32 = -t128 * t185 + t170;
t165 = t128 * t181;
t34 = t126 * t185 + t165;
t22 = t131 * t32 - t133 * t34;
t225 = pkin(6) * t22;
t132 = sin(qJ(2));
t134 = cos(qJ(2));
t158 = t96 * qJD(3);
t152 = qJD(2) * qJD(3);
t148 = t133 * t152;
t151 = t131 * qJDD(2);
t103 = t148 + t151;
t149 = t131 * t152;
t150 = t133 * qJDD(2);
t142 = -t149 + t150;
t73 = t128 * t103 + t126 * t142;
t198 = t158 - t73;
t224 = t132 * t22 + t134 * t198;
t223 = pkin(3) * t32;
t222 = qJ(4) * t32;
t221 = qJ(4) * t34;
t180 = t96 ^ 2;
t82 = t180 - t135;
t220 = t131 * (-t128 * t82 + t170) - t133 * (t126 * t82 + t165);
t219 = t198 * qJ(5);
t182 = qJDD(3) - t175;
t169 = t126 * t182;
t183 = -t180 - t135;
t188 = t128 * t183 - t169;
t57 = t128 * t182;
t191 = t126 * t183 + t57;
t201 = -t131 * t191 + t133 * t188;
t145 = t126 * t103 - t128 * t142;
t157 = t98 * qJD(3);
t48 = t145 + t157;
t216 = -pkin(2) * t48 + pkin(6) * t201;
t50 = -t145 + t157;
t53 = t158 + t73;
t190 = t126 * t53 + t128 * t50;
t192 = t126 * t50 - t128 * t53;
t200 = -t131 * t192 + t133 * t190;
t47 = t95 + t180;
t215 = pkin(2) * t47 + pkin(6) * t200;
t214 = t132 * t201 - t134 * t48;
t213 = t132 * t200 + t134 * t47;
t211 = pkin(3) * t192;
t207 = qJ(4) * t188;
t206 = qJ(4) * t191;
t205 = qJ(4) * t192;
t162 = qJD(4) * t98;
t204 = pkin(3) * t191 - 0.2e1 * t162;
t203 = pkin(3) * t47 + qJ(4) * t190;
t202 = t131 * (-t126 * t198 + t128 * t48) - t133 * (-t126 * t48 - t128 * t198);
t83 = -t95 + t135;
t199 = t131 * (-t126 * t83 + t57) + t133 * (t128 * t83 + t169);
t184 = t95 - t180;
t136 = qJD(2) ^ 2;
t113 = t131 * t136 * t133;
t110 = qJDD(3) + t113;
t127 = sin(pkin(7));
t129 = cos(pkin(7));
t107 = -t127 * g(1) + t129 * g(2);
t108 = -t129 * g(1) - t127 * g(2);
t156 = -g(3) + qJDD(1);
t81 = t134 * t108 + t132 * t156;
t75 = -t136 * pkin(2) + qJDD(2) * pkin(6) + t81;
t55 = -t133 * t107 + t131 * t75;
t137 = (-t103 + t148) * qJ(4) + t110 * pkin(3) - t55;
t109 = qJD(3) * pkin(3) - qJ(4) * t155;
t123 = t133 ^ 2;
t118 = t123 * t136;
t56 = t131 * t107 + t133 * t75;
t36 = -pkin(3) * t118 + t142 * qJ(4) - qJD(3) * t109 + t56;
t18 = -0.2e1 * qJD(4) * t96 + t126 * t137 + t128 * t36;
t179 = 2 * qJD(5);
t178 = pkin(4) * t126;
t177 = pkin(4) * t128;
t146 = t126 * t36 - t128 * t137;
t17 = t146 + 0.2e1 * t162;
t5 = t126 * t18 - t128 * t17;
t176 = t131 * t5;
t80 = -t132 * t108 + t134 * t156;
t74 = -qJDD(2) * pkin(2) - t136 * pkin(6) - t80;
t41 = -t142 * pkin(3) - qJ(4) * t118 + t109 * t155 + qJDD(4) + t74;
t173 = t126 * t41;
t167 = t128 * t41;
t164 = qJ(5) * t128;
t161 = t131 * t110;
t159 = t133 * (qJDD(3) - t113);
t154 = qJD(3) * t126;
t153 = qJD(3) * t128;
t147 = -qJ(5) * t126 - pkin(3);
t6 = t126 * t17 + t128 * t18;
t28 = t131 * t55 + t133 * t56;
t143 = -qJDD(3) * pkin(4) - t135 * qJ(5) + qJDD(5) + t146;
t58 = t96 * pkin(4) - t98 * qJ(5);
t10 = (0.2e1 * qJD(4) + t58) * t98 + t143;
t144 = qJDD(3) * qJ(5) + qJD(3) * t179 - t96 * t58 + t18;
t9 = -t135 * pkin(4) + t144;
t3 = -t128 * t10 + t126 * t9;
t1 = -t131 * t3 + t133 * (t126 * t10 + t128 * t9);
t104 = -0.2e1 * t149 + t150;
t141 = t131 * (t126 * t145 + t96 * t153) + t133 * (-t128 * t145 + t96 * t154);
t79 = t98 * t154;
t140 = t131 * t79 + (-t96 * t160 + t133 * (-t126 * t96 - t128 * t98)) * qJD(3);
t139 = t145 * pkin(4) + t219 + t41;
t138 = t98 * t179 - t139;
t122 = t131 ^ 2;
t117 = t122 * t136;
t106 = t117 + t118;
t105 = (t122 + t123) * qJDD(2);
t102 = 0.2e1 * t148 + t151;
t77 = -t159 - t131 * (-t117 - t135);
t76 = t133 * (-t118 - t135) - t161;
t23 = t131 * (t128 * t73 - t79) + t133 * (t126 * t73 + t98 * t153);
t15 = (pkin(4) * qJD(3) - (2 * qJD(5))) * t98 + t139;
t12 = (-t48 - t157) * pkin(4) + t138;
t11 = -pkin(4) * t157 + t138 - t219;
t8 = qJ(5) * t47 + t10;
t7 = (-t135 + t47) * pkin(4) + t144;
t2 = t133 * t6 - t176;
t4 = [0, 0, 0, 0, 0, 0, 0, 0, 0, t156, 0, 0, 0, 0, 0, 0, t134 * qJDD(2) - t132 * t136, -t132 * qJDD(2) - t134 * t136, 0, t132 * t81 + t134 * t80, 0, 0, 0, 0, 0, 0, t134 * t104 + t132 * t76, -t134 * t102 + t132 * t77, t132 * t105 + t134 * t106, t132 * t28 - t134 * t74, 0, 0, 0, 0, 0, 0, t214, t224, t213, t132 * t2 - t134 * t41, 0, 0, 0, 0, 0, 0, t214, t213, -t224, t132 * t1 - t134 * t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(2), t80, -t81, 0, 0, (t103 + t148) * t131, t133 * t102 + t131 * t104, t161 + t133 * (-t117 + t135), t104 * t133, t131 * (t118 - t135) + t159, 0, pkin(2) * t104 + pkin(6) * t76 - t133 * t74, -pkin(2) * t102 + pkin(6) * t77 + t131 * t74, pkin(2) * t106 + pkin(6) * t105 + t28, -pkin(2) * t74 + pkin(6) * t28, t23, -t202, t199, t141, -t220, t140, t131 * (t173 - t206) + t133 * (-pkin(3) * t48 - t167 + t207) + t216, t131 * (t167 + t222) + t133 * (pkin(3) * t198 + t173 - t221) + pkin(2) * t198 + t225, t131 * (-t5 - t205) + t133 * (t203 + t6) + t215, -qJ(4) * t176 + t133 * (-pkin(3) * t41 + qJ(4) * t6) - pkin(2) * t41 + pkin(6) * t2, t23, t199, t202, t140, t220, t141, t131 * (-t126 * t12 - t48 * t164 - t206) + t133 * (t128 * t12 + t147 * t48 + t207) + t216, t131 * (-t126 * t7 + t128 * t8 - t205) + t133 * (t126 * t8 + t128 * t7 + t203) + t215, t131 * (t128 * t11 - t222) + t133 * (t126 * t11 + t221) - t225 - (-t131 * t178 + t133 * (pkin(3) + t177) + pkin(2)) * t198, (t131 * (-t164 + t178) + t133 * (t147 - t177) - pkin(2)) * t15 + (pkin(6) + qJ(4)) * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t113, t117 - t118, t151, t113, t150, qJDD(3), -t55, -t56, 0, 0, t175, t184, t53, -t175, t50, qJDD(3), -t146 + t204, -t18 - t223, t211, pkin(3) * t5, t175, t53, -t184, qJDD(3), -t50, -t175, pkin(4) * t182 + qJ(5) * t183 - t98 * t58 - t143 + t204, -pkin(4) * t53 + qJ(5) * t50 + t211, t223 + qJ(5) * t181 + (-t135 - t185) * pkin(4) + t144, pkin(3) * t3 - pkin(4) * t10 + qJ(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t48, -t198, -t47, t41, 0, 0, 0, 0, 0, 0, t48, -t47, t198, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t182, t53, t185, t10;];
tauJ_reg = t4;
