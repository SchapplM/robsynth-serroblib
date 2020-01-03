% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRPPR1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2020-01-03 11:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRPPR1_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_invdynJ_fixb_reg2_snew_vp: pkin has to be [9x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:56:00
% EndTime: 2020-01-03 11:56:04
% DurationCPUTime: 1.20s
% Computational Cost: add. (5429->175), mult. (7650->271), div. (0->0), fcn. (4733->10), ass. (0->124)
t131 = (qJD(1) + qJD(2));
t129 = t131 ^ 2;
t130 = qJDD(1) + qJDD(2);
t134 = sin(pkin(8));
t136 = cos(pkin(8));
t139 = sin(qJ(1));
t142 = cos(qJ(1));
t155 = -t142 * g(2) - t139 * g(3);
t112 = qJDD(1) * pkin(1) + t155;
t154 = t139 * g(2) - t142 * g(3);
t113 = -qJD(1) ^ 2 * pkin(1) - t154;
t138 = sin(qJ(2));
t141 = cos(qJ(2));
t87 = t141 * t112 - t138 * t113;
t149 = t130 * pkin(2) + t87;
t88 = t138 * t112 + t141 * t113;
t82 = -t129 * pkin(2) + t88;
t60 = t134 * t149 + t136 * t82;
t190 = -t129 * pkin(3) + t130 * qJ(4) + (2 * qJD(4) * t131) + t60;
t137 = sin(qJ(5));
t133 = sin(pkin(9));
t135 = cos(pkin(9));
t140 = cos(qJ(5));
t187 = -t133 * t137 + t135 * t140;
t97 = t187 * t131;
t152 = t133 * t140 + t135 * t137;
t99 = t152 * t131;
t76 = t99 * t97;
t183 = qJDD(5) + t76;
t189 = t137 * t183;
t188 = t140 * t183;
t132 = -g(1) + qJDD(3);
t126 = t135 * t132;
t46 = t190 * t133 - t126;
t47 = t133 * t132 + t190 * t135;
t22 = t133 * t46 + t135 * t47;
t186 = t129 * t135;
t144 = t133 ^ 2;
t146 = t135 ^ 2;
t185 = t144 + t146;
t184 = t126 + (pkin(4) * t186 - pkin(7) * t130 - t190) * t133;
t111 = t185 * t129;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t171 = t146 * t129;
t174 = t135 * t130;
t41 = -pkin(4) * t171 + pkin(7) * t174 + t47;
t14 = t137 * t41 - t140 * t184;
t15 = t184 * t137 + t140 * t41;
t8 = t137 * t15 - t140 * t14;
t182 = t133 * t8;
t59 = -t134 * t82 + t136 * t149;
t55 = -t130 * pkin(3) - t129 * qJ(4) + qJDD(4) - t59;
t43 = -pkin(4) * t174 + t55 + (-t129 * t144 - t171) * pkin(7);
t181 = t137 * t43;
t70 = qJDD(5) - t76;
t180 = t137 * t70;
t179 = t140 * t43;
t178 = t140 * t70;
t177 = t133 * t130;
t175 = t134 * t130;
t172 = t136 * t130;
t170 = t97 * qJD(5);
t169 = t99 * qJD(5);
t166 = qJD(5) * t137;
t165 = qJD(5) * t140;
t11 = t134 * t22 - t136 * t55;
t164 = pkin(2) * t11 - pkin(3) * t55 + qJ(4) * t22;
t109 = -t136 * t129 - t175;
t162 = pkin(2) * t109 - t60;
t9 = t137 * t14 + t140 * t15;
t103 = t185 * t186;
t86 = -t134 * t103 + t135 * t172;
t161 = pkin(2) * t86 + pkin(3) * t174 - qJ(4) * t103 - t135 * t55;
t63 = t187 * t130;
t94 = t152 * t130;
t56 = t137 * t63 - t140 * t94;
t57 = t137 * t94 + t140 * t63;
t31 = -t133 * t56 + t135 * t57;
t65 = -t95 - t96;
t17 = t134 * t31 - t136 * t65;
t160 = t133 * (-pkin(7) * t56 - t8) + t135 * (-pkin(4) * t65 + pkin(7) * t57 + t9) - pkin(3) * t65 + qJ(4) * t31 + pkin(2) * t17;
t143 = qJD(5) ^ 2;
t68 = -t143 - t95;
t53 = t137 * t68 + t188;
t54 = t140 * t68 - t189;
t29 = -t133 * t53 + t135 * t54;
t72 = -t63 + 0.2e1 * t169;
t19 = t134 * t29 - t136 * t72;
t159 = t133 * (-pkin(7) * t53 + t181) + t135 * (-pkin(4) * t72 + pkin(7) * t54 - t179) - pkin(3) * t72 + qJ(4) * t29 + pkin(2) * t19;
t91 = -t96 - t143;
t61 = t140 * t91 - t180;
t62 = -t137 * t91 - t178;
t37 = -t133 * t61 + t135 * t62;
t74 = t94 + 0.2e1 * t170;
t25 = t134 * t37 - t136 * t74;
t158 = t133 * (-pkin(7) * t61 + t179) + t135 * (-pkin(4) * t74 + pkin(7) * t62 + t181) - pkin(3) * t74 + qJ(4) * t37 + pkin(2) * t25;
t153 = t134 * t129 - t172;
t157 = -pkin(2) * t153 + t59;
t123 = t144 * t130;
t124 = t146 * t130;
t108 = t124 + t123;
t81 = t134 * t108 + t136 * t111;
t156 = pkin(2) * t81 + pkin(3) * t111 + qJ(4) * t108 + t22;
t102 = t133 * t111;
t85 = t134 * t102 - t133 * t172;
t151 = pkin(2) * t85 - pkin(3) * t177 + qJ(4) * t102 + t133 * t55;
t4 = t135 * t9 - t182;
t2 = t134 * t4 - t136 * t43;
t150 = pkin(2) * t2 - pkin(7) * t182 + qJ(4) * t4 - pkin(3) * t43 + t135 * (-pkin(4) * t43 + pkin(7) * t9);
t114 = 0.2e1 * t133 * t174;
t90 = -t96 + t143;
t89 = t95 - t143;
t75 = t94 + t170;
t73 = t63 - t169;
t48 = (t133 * (t137 * t99 + t140 * t97) + t135 * (t137 * t97 - t140 * t99)) * qJD(5);
t39 = t133 * (t140 * t75 - t99 * t166) + t135 * (t137 * t75 + t99 * t165);
t38 = t133 * (-t137 * t73 - t97 * t165) + t135 * (t140 * t73 - t97 * t166);
t36 = t133 * (-t137 * t90 + t188) + t135 * (t140 * t90 + t189);
t35 = t133 * (t140 * t89 - t180) + t135 * (t137 * t89 + t178);
t33 = t134 * t60 + t136 * t59;
t32 = pkin(2) * t33;
t30 = t133 * (-t137 * t74 - t140 * t72) + t135 * (-t137 * t72 + t140 * t74);
t1 = [0, 0, 0, 0, 0, qJDD(1), t155, t154, 0, 0, 0, 0, 0, 0, 0, t130, pkin(1) * (-t138 * t129 + t141 * t130) + t87, pkin(1) * (-t141 * t129 - t138 * t130) - t88, 0, pkin(1) * (t138 * t88 + t141 * t87), 0, 0, 0, 0, 0, t130, pkin(1) * (t138 * t109 - t141 * t153) + t157, pkin(1) * (t141 * t109 + t138 * t153) + t162, 0, pkin(1) * (t138 * (-t134 * t59 + t136 * t60) + t141 * t33) + t32, t123, t114, 0, t124, 0, 0, pkin(1) * (t138 * (-t136 * t103 - t134 * t174) + t141 * t86) + t161, pkin(1) * (t138 * (t136 * t102 + t133 * t175) + t141 * t85) + t151, pkin(1) * (t138 * (t136 * t108 - t134 * t111) + t141 * t81) + t156, pkin(1) * (t138 * (t134 * t55 + t136 * t22) + t141 * t11) + t164, t39, t30, t36, t38, t35, t48, pkin(1) * (t138 * (t134 * t72 + t136 * t29) + t141 * t19) + t159, pkin(1) * (t138 * (t134 * t74 + t136 * t37) + t141 * t25) + t158, pkin(1) * (t138 * (t134 * t65 + t136 * t31) + t141 * t17) + t160, pkin(1) * (t138 * (t134 * t43 + t136 * t4) + t141 * t2) + t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t130, t87, -t88, 0, 0, 0, 0, 0, 0, 0, t130, t157, t162, 0, t32, t123, t114, 0, t124, 0, 0, t161, t151, t156, t164, t39, t30, t36, t38, t35, t48, t159, t158, t160, t150; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133 * t47 - t135 * t46, 0, 0, 0, 0, 0, 0, t133 * t54 + t135 * t53, t133 * t62 + t135 * t61, t133 * t57 + t135 * t56, t133 * t9 + t135 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t174, t177, -t111, t55, 0, 0, 0, 0, 0, 0, t72, t74, t65, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t96 - t95, t94, t76, t63, qJDD(5), -t14, -t15, 0, 0;];
tauJ_reg = t1;
