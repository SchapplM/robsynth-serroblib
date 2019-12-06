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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:18:28
% EndTime: 2019-12-05 18:18:35
% DurationCPUTime: 1.25s
% Computational Cost: add. (5429->175), mult. (7650->271), div. (0->0), fcn. (4733->10), ass. (0->124)
t133 = (qJD(1) + qJD(2));
t131 = t133 ^ 2;
t132 = qJDD(1) + qJDD(2);
t136 = sin(pkin(8));
t138 = cos(pkin(8));
t141 = sin(qJ(1));
t144 = cos(qJ(1));
t169 = t144 * g(2) + t141 * g(3);
t112 = qJDD(1) * pkin(1) + t169;
t156 = -t141 * g(2) + t144 * g(3);
t113 = -qJD(1) ^ 2 * pkin(1) - t156;
t140 = sin(qJ(2));
t143 = cos(qJ(2));
t87 = t143 * t112 - t140 * t113;
t151 = t132 * pkin(2) + t87;
t88 = t140 * t112 + t143 * t113;
t82 = -t131 * pkin(2) + t88;
t60 = t136 * t151 + t138 * t82;
t192 = -t131 * pkin(3) + t132 * qJ(4) + (2 * qJD(4) * t133) + t60;
t139 = sin(qJ(5));
t135 = sin(pkin(9));
t137 = cos(pkin(9));
t142 = cos(qJ(5));
t189 = -t135 * t139 + t137 * t142;
t97 = t189 * t133;
t154 = t135 * t142 + t137 * t139;
t99 = t154 * t133;
t76 = t99 * t97;
t185 = qJDD(5) + t76;
t191 = t139 * t185;
t190 = t142 * t185;
t134 = -g(1) + qJDD(3);
t126 = t137 * t134;
t46 = t192 * t135 - t126;
t47 = t135 * t134 + t192 * t137;
t22 = t135 * t46 + t137 * t47;
t188 = t131 * t137;
t146 = t135 ^ 2;
t148 = t137 ^ 2;
t187 = t146 + t148;
t186 = t126 + (pkin(4) * t188 - pkin(7) * t132 - t192) * t135;
t111 = t187 * t131;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t173 = t148 * t131;
t176 = t137 * t132;
t41 = -pkin(4) * t173 + pkin(7) * t176 + t47;
t14 = t139 * t41 - t142 * t186;
t15 = t186 * t139 + t142 * t41;
t8 = t139 * t15 - t142 * t14;
t184 = t135 * t8;
t59 = -t136 * t82 + t138 * t151;
t55 = -t132 * pkin(3) - t131 * qJ(4) + qJDD(4) - t59;
t43 = -pkin(4) * t176 + t55 + (-t131 * t146 - t173) * pkin(7);
t183 = t139 * t43;
t70 = qJDD(5) - t76;
t182 = t139 * t70;
t181 = t142 * t43;
t180 = t142 * t70;
t179 = t135 * t132;
t177 = t136 * t132;
t174 = t138 * t132;
t172 = t97 * qJD(5);
t171 = t99 * qJD(5);
t167 = qJD(5) * t139;
t166 = qJD(5) * t142;
t11 = t136 * t22 - t138 * t55;
t165 = pkin(2) * t11 - pkin(3) * t55 + qJ(4) * t22;
t109 = -t138 * t131 - t177;
t163 = pkin(2) * t109 - t60;
t9 = t139 * t14 + t142 * t15;
t103 = t187 * t188;
t86 = -t136 * t103 + t137 * t174;
t162 = pkin(2) * t86 + pkin(3) * t176 - qJ(4) * t103 - t137 * t55;
t63 = t189 * t132;
t94 = t154 * t132;
t56 = t139 * t63 - t142 * t94;
t57 = t139 * t94 + t142 * t63;
t31 = -t135 * t56 + t137 * t57;
t65 = -t95 - t96;
t17 = t136 * t31 - t138 * t65;
t161 = t135 * (-pkin(7) * t56 - t8) + t137 * (-pkin(4) * t65 + pkin(7) * t57 + t9) - pkin(3) * t65 + qJ(4) * t31 + pkin(2) * t17;
t145 = qJD(5) ^ 2;
t68 = -t145 - t95;
t53 = t139 * t68 + t190;
t54 = t142 * t68 - t191;
t29 = -t135 * t53 + t137 * t54;
t72 = -t63 + 0.2e1 * t171;
t19 = t136 * t29 - t138 * t72;
t160 = t135 * (-pkin(7) * t53 + t183) + t137 * (-pkin(4) * t72 + pkin(7) * t54 - t181) - pkin(3) * t72 + qJ(4) * t29 + pkin(2) * t19;
t91 = -t96 - t145;
t61 = t142 * t91 - t182;
t62 = -t139 * t91 - t180;
t37 = -t135 * t61 + t137 * t62;
t74 = t94 + 0.2e1 * t172;
t25 = t136 * t37 - t138 * t74;
t159 = t135 * (-pkin(7) * t61 + t181) + t137 * (-pkin(4) * t74 + pkin(7) * t62 + t183) - pkin(3) * t74 + qJ(4) * t37 + pkin(2) * t25;
t155 = t136 * t131 - t174;
t158 = -pkin(2) * t155 + t59;
t123 = t146 * t132;
t124 = t148 * t132;
t108 = t124 + t123;
t81 = t136 * t108 + t138 * t111;
t157 = pkin(2) * t81 + pkin(3) * t111 + qJ(4) * t108 + t22;
t102 = t135 * t111;
t85 = t136 * t102 - t135 * t174;
t153 = pkin(2) * t85 - pkin(3) * t179 + qJ(4) * t102 + t135 * t55;
t4 = t137 * t9 - t184;
t2 = t136 * t4 - t138 * t43;
t152 = pkin(2) * t2 - pkin(7) * t184 + qJ(4) * t4 - pkin(3) * t43 + t137 * (-pkin(4) * t43 + pkin(7) * t9);
t114 = 0.2e1 * t135 * t176;
t90 = -t96 + t145;
t89 = t95 - t145;
t75 = t94 + t172;
t73 = t63 - t171;
t48 = (t135 * (t139 * t99 + t142 * t97) + t137 * (t139 * t97 - t142 * t99)) * qJD(5);
t39 = t135 * (t142 * t75 - t99 * t167) + t137 * (t139 * t75 + t99 * t166);
t38 = t135 * (-t139 * t73 - t97 * t166) + t137 * (t142 * t73 - t97 * t167);
t36 = t135 * (-t139 * t90 + t190) + t137 * (t142 * t90 + t191);
t35 = t135 * (t142 * t89 - t182) + t137 * (t139 * t89 + t180);
t33 = t136 * t60 + t138 * t59;
t32 = pkin(2) * t33;
t30 = t135 * (-t139 * t74 - t142 * t72) + t137 * (-t139 * t72 + t142 * t74);
t1 = [0, 0, 0, 0, 0, qJDD(1), t169, t156, 0, 0, 0, 0, 0, 0, 0, t132, pkin(1) * (-t140 * t131 + t143 * t132) + t87, pkin(1) * (-t143 * t131 - t140 * t132) - t88, 0, pkin(1) * (t140 * t88 + t143 * t87), 0, 0, 0, 0, 0, t132, pkin(1) * (t140 * t109 - t143 * t155) + t158, pkin(1) * (t143 * t109 + t140 * t155) + t163, 0, pkin(1) * (t140 * (-t136 * t59 + t138 * t60) + t143 * t33) + t32, t123, t114, 0, t124, 0, 0, pkin(1) * (t140 * (-t138 * t103 - t136 * t176) + t143 * t86) + t162, pkin(1) * (t140 * (t138 * t102 + t135 * t177) + t143 * t85) + t153, pkin(1) * (t140 * (t138 * t108 - t136 * t111) + t143 * t81) + t157, pkin(1) * (t140 * (t136 * t55 + t138 * t22) + t143 * t11) + t165, t39, t30, t36, t38, t35, t48, pkin(1) * (t140 * (t136 * t72 + t138 * t29) + t143 * t19) + t160, pkin(1) * (t140 * (t136 * t74 + t138 * t37) + t143 * t25) + t159, pkin(1) * (t140 * (t136 * t65 + t138 * t31) + t143 * t17) + t161, pkin(1) * (t140 * (t136 * t43 + t138 * t4) + t143 * t2) + t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t132, t87, -t88, 0, 0, 0, 0, 0, 0, 0, t132, t158, t163, 0, t32, t123, t114, 0, t124, 0, 0, t162, t153, t157, t165, t39, t30, t36, t38, t35, t48, t160, t159, t161, t152; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134, 0, 0, 0, 0, 0, 0, 0, 0, 0, t135 * t47 - t137 * t46, 0, 0, 0, 0, 0, 0, t135 * t54 + t137 * t53, t135 * t62 + t137 * t61, t135 * t57 + t137 * t56, t135 * t9 + t137 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t176, t179, -t111, t55, 0, 0, 0, 0, 0, 0, t72, t74, t65, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t96 - t95, t94, t76, t63, qJDD(5), -t14, -t15, 0, 0;];
tauJ_reg = t1;
