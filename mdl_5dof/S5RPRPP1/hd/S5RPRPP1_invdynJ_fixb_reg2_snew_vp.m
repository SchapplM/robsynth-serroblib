% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRPP1
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
%   pkin=[a2,a3,a4,a5,d1,d3,theta2,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:09
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRPP1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPP1_invdynJ_fixb_reg2_snew_vp: pkin has to be [8x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:09:14
% EndTime: 2019-12-31 18:09:19
% DurationCPUTime: 1.88s
% Computational Cost: add. (3389->235), mult. (7447->297), div. (0->0), fcn. (4667->8), ass. (0->142)
t130 = cos(pkin(8));
t128 = sin(pkin(8));
t134 = cos(qJ(3));
t133 = sin(qJ(3));
t162 = qJD(1) * t133;
t94 = -t130 * t134 * qJD(1) + t128 * t162;
t165 = t130 * t133;
t96 = (t128 * t134 + t165) * qJD(1);
t182 = t96 * t94;
t190 = qJDD(3) + t182;
t176 = t128 * t190;
t135 = qJD(3) ^ 2;
t93 = t96 ^ 2;
t194 = -t93 - t135;
t34 = -t130 * t194 + t176;
t232 = pkin(3) * t34;
t231 = qJ(4) * t34;
t171 = t130 * t190;
t36 = t128 * t194 + t171;
t230 = qJ(4) * t36;
t229 = t133 * t36 + t134 * t34;
t22 = t133 * t34 - t134 * t36;
t189 = t94 ^ 2;
t80 = t189 - t135;
t228 = t133 * (-t130 * t80 + t176) - t134 * (t128 * t80 + t171);
t169 = qJD(3) * t94;
t159 = qJD(1) * qJD(3);
t153 = t134 * t159;
t158 = t133 * qJDD(1);
t104 = t153 + t158;
t154 = t133 * t159;
t157 = t134 * qJDD(1);
t145 = -t154 + t157;
t75 = t130 * t104 + t128 * t145;
t206 = t169 - t75;
t227 = t206 * qJ(5);
t129 = sin(pkin(7));
t131 = cos(pkin(7));
t191 = qJDD(3) - t182;
t175 = t128 * t191;
t192 = -t189 - t135;
t197 = t130 * t192 - t175;
t57 = t130 * t191;
t200 = t128 * t192 + t57;
t211 = -t133 * t200 + t134 * t197;
t149 = t104 * t128 - t130 * t145;
t168 = qJD(3) * t96;
t49 = t149 + t168;
t224 = pkin(1) * (t129 * t211 - t131 * t49) + pkin(6) * t211 - pkin(2) * t49;
t51 = -t149 + t168;
t54 = t169 + t75;
t199 = t128 * t54 + t130 * t51;
t201 = t128 * t51 - t130 * t54;
t209 = -t133 * t201 + t134 * t199;
t48 = t93 + t189;
t223 = pkin(2) * t48 + pkin(1) * (t129 * t209 + t131 * t48) + pkin(6) * t209;
t221 = pkin(3) * t201;
t217 = qJ(4) * t197;
t216 = qJ(4) * t200;
t215 = qJ(4) * t201;
t166 = qJD(4) * t96;
t214 = pkin(3) * t200 - 0.2e1 * t166;
t213 = pkin(3) * t48 + qJ(4) * t199;
t212 = t133 * (-t128 * t206 + t130 * t49) - t134 * (-t128 * t49 - t130 * t206);
t210 = t133 * t197 + t134 * t200;
t208 = t133 * t199 + t134 * t201;
t81 = -t93 + t135;
t207 = t133 * (-t128 * t81 + t57) + t134 * (t130 * t81 + t175);
t193 = t93 - t189;
t136 = qJD(1) ^ 2;
t114 = t134 * t136 * t133;
t109 = qJDD(3) + t114;
t125 = -g(3) + qJDD(2);
t186 = sin(qJ(1));
t187 = cos(qJ(1));
t144 = t187 * g(1) + t186 * g(2);
t102 = -t136 * pkin(1) - t144;
t143 = t186 * g(1) - t187 * g(2);
t142 = qJDD(1) * pkin(1) + t143;
t180 = t131 * t102 + t129 * t142;
t59 = -pkin(2) * t136 + qJDD(1) * pkin(6) + t180;
t55 = -t134 * t125 + t133 * t59;
t137 = (-t104 + t153) * qJ(4) + t109 * pkin(3) - t55;
t108 = qJD(3) * pkin(3) - qJ(4) * t162;
t124 = t134 ^ 2;
t119 = t124 * t136;
t56 = t133 * t125 + t134 * t59;
t32 = -pkin(3) * t119 + t145 * qJ(4) - qJD(3) * t108 + t56;
t18 = -0.2e1 * qJD(4) * t94 + t128 * t137 + t130 * t32;
t188 = 2 * qJD(5);
t185 = pkin(4) * t128;
t184 = pkin(4) * t130;
t151 = t128 * t32 - t130 * t137;
t17 = t151 + 0.2e1 * t166;
t5 = t128 * t18 - t130 * t17;
t183 = t133 * t5;
t150 = -t129 * t102 + t131 * t142;
t58 = -qJDD(1) * pkin(2) - t136 * pkin(6) - t150;
t37 = -t145 * pkin(3) - qJ(4) * t119 + t108 * t162 + qJDD(4) + t58;
t179 = t128 * t37;
t173 = t130 * t37;
t170 = qJ(5) * t130;
t164 = t133 * t109;
t110 = qJDD(3) - t114;
t163 = t134 * t110;
t161 = qJD(3) * t128;
t160 = qJD(3) * t130;
t156 = -pkin(1) * t131 - pkin(2);
t155 = pkin(1) * t129 + pkin(6);
t152 = -qJ(5) * t128 - pkin(3);
t6 = t128 * t17 + t130 * t18;
t28 = t133 * t55 + t134 * t56;
t60 = pkin(4) * t94 - qJ(5) * t96;
t147 = qJDD(3) * qJ(5) + qJD(3) * t188 - t94 * t60 + t18;
t146 = -qJDD(3) * pkin(4) - t135 * qJ(5) + qJDD(5) + t151;
t105 = -0.2e1 * t154 + t157;
t10 = (0.2e1 * qJD(4) + t60) * t96 + t146;
t141 = t133 * (t128 * t149 + t94 * t160) + t134 * (-t130 * t149 + t94 * t161);
t140 = t149 * pkin(4) + t227 + t37;
t79 = t96 * t161;
t139 = t133 * t79 + (-t94 * t165 + t134 * (-t128 * t94 - t130 * t96)) * qJD(3);
t138 = t96 * t188 - t140;
t123 = t133 ^ 2;
t118 = t123 * t136;
t113 = -t119 - t135;
t112 = -t118 - t135;
t107 = t118 + t119;
t106 = (t123 + t124) * qJDD(1);
t103 = 0.2e1 * t153 + t158;
t77 = -t112 * t133 - t163;
t76 = t113 * t134 - t164;
t23 = t133 * (t130 * t75 - t79) + t134 * (t128 * t75 + t96 * t160);
t15 = (pkin(4) * qJD(3) - (2 * qJD(5))) * t96 + t140;
t12 = (-t49 - t168) * pkin(4) + t138;
t11 = -pkin(4) * t168 + t138 - t227;
t9 = -pkin(4) * t135 + t147;
t8 = qJ(5) * t48 + t10;
t7 = (-t135 + t48) * pkin(4) + t147;
t4 = t10 * t128 + t130 * t9;
t3 = -t10 * t130 + t128 * t9;
t2 = t134 * t6 - t183;
t1 = [0, 0, 0, 0, 0, qJDD(1), t143, t144, 0, 0, 0, 0, 0, 0, 0, qJDD(1), pkin(1) * (qJDD(1) * t131 - t129 * t136) + t150, pkin(1) * (-qJDD(1) * t129 - t131 * t136) - t180, 0, pkin(1) * (t129 * t180 + t131 * t150), (t104 + t153) * t133, t103 * t134 + t105 * t133, t164 + t134 * (-t118 + t135), t105 * t134, t133 * (t119 - t135) + t163, 0, -t134 * t58 + pkin(2) * t105 + pkin(6) * t76 + pkin(1) * (t105 * t131 + t129 * t76), t133 * t58 - pkin(2) * t103 + pkin(6) * t77 + pkin(1) * (-t103 * t131 + t129 * t77), pkin(2) * t107 + pkin(6) * t106 + pkin(1) * (t106 * t129 + t107 * t131) + t28, -pkin(2) * t58 + pkin(6) * t28 + pkin(1) * (t129 * t28 - t131 * t58), t23, -t212, t207, t141, -t228, t139, t133 * (t179 - t216) + t134 * (-pkin(3) * t49 - t173 + t217) + t224, t133 * (t173 + t231) + t134 * (pkin(3) * t206 + t179 - t230) + pkin(2) * t206 + pkin(6) * t22 + pkin(1) * (t129 * t22 + t131 * t206), t133 * (-t5 - t215) + t134 * (t213 + t6) + t223, -qJ(4) * t183 + t134 * (-pkin(3) * t37 + qJ(4) * t6) - pkin(2) * t37 + pkin(6) * t2 + pkin(1) * (t129 * t2 - t131 * t37), t23, t207, t212, t139, t228, t141, t133 * (-t12 * t128 - t49 * t170 - t216) + t134 * (t12 * t130 + t152 * t49 + t217) + t224, t133 * (-t128 * t7 + t130 * t8 - t215) + t134 * (t128 * t8 + t130 * t7 + t213) + t223, t133 * (t11 * t130 - t231) + t134 * (t11 * t128 + t230) - t155 * t22 - (-t133 * t185 + t134 * (pkin(3) + t184) - t156) * t206, (t133 * (-t170 + t185) + t134 * (t152 - t184) + t156) * t15 + (t155 + qJ(4)) * (-t133 * t3 + t134 * t4); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t125, 0, 0, 0, 0, 0, 0, t109 * t134 + t113 * t133, -t110 * t133 + t112 * t134, 0, t133 * t56 - t134 * t55, 0, 0, 0, 0, 0, 0, t210, -t229, t208, t133 * t6 + t134 * t5, 0, 0, 0, 0, 0, 0, t210, t208, t229, t133 * t4 + t134 * t3; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t114, t118 - t119, t158, t114, t157, qJDD(3), -t55, -t56, 0, 0, t182, t193, t54, -t182, t51, qJDD(3), -t151 + t214, -t18 - t232, t221, pkin(3) * t5, t182, t54, -t193, qJDD(3), -t51, -t182, pkin(4) * t191 + qJ(5) * t192 - t60 * t96 - t146 + t214, -pkin(4) * t54 + qJ(5) * t51 + t221, t232 + qJ(5) * t190 + (-t135 - t194) * pkin(4) + t147, pkin(3) * t3 - pkin(4) * t10 + qJ(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t49, -t206, -t48, t37, 0, 0, 0, 0, 0, 0, t49, -t48, t206, t15; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t191, t54, t194, t10;];
tauJ_reg = t1;
