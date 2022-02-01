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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:43
% EndTime: 2022-01-20 09:51:48
% DurationCPUTime: 1.27s
% Computational Cost: add. (5429->175), mult. (7650->271), div. (0->0), fcn. (4733->10), ass. (0->124)
t132 = (qJD(1) + qJD(2));
t130 = t132 ^ 2;
t131 = qJDD(1) + qJDD(2);
t135 = sin(pkin(8));
t137 = cos(pkin(8));
t140 = sin(qJ(1));
t143 = cos(qJ(1));
t162 = t140 * g(1) - t143 * g(2);
t112 = qJDD(1) * pkin(1) + t162;
t155 = t143 * g(1) + t140 * g(2);
t113 = -qJD(1) ^ 2 * pkin(1) - t155;
t139 = sin(qJ(2));
t142 = cos(qJ(2));
t87 = t142 * t112 - t139 * t113;
t150 = t131 * pkin(2) + t87;
t88 = t139 * t112 + t142 * t113;
t82 = -t130 * pkin(2) + t88;
t60 = t135 * t150 + t137 * t82;
t191 = -t130 * pkin(3) + t131 * qJ(4) + (2 * qJD(4) * t132) + t60;
t138 = sin(qJ(5));
t134 = sin(pkin(9));
t136 = cos(pkin(9));
t141 = cos(qJ(5));
t188 = -t134 * t138 + t136 * t141;
t97 = t188 * t132;
t153 = t134 * t141 + t136 * t138;
t99 = t153 * t132;
t76 = t99 * t97;
t184 = qJDD(5) + t76;
t190 = t138 * t184;
t189 = t141 * t184;
t133 = -g(3) + qJDD(3);
t126 = t136 * t133;
t46 = t191 * t134 - t126;
t47 = t134 * t133 + t191 * t136;
t22 = t134 * t46 + t136 * t47;
t187 = t130 * t136;
t145 = t134 ^ 2;
t147 = t136 ^ 2;
t186 = t145 + t147;
t185 = t126 + (pkin(4) * t187 - pkin(7) * t131 - t191) * t134;
t111 = t186 * t130;
t95 = t97 ^ 2;
t96 = t99 ^ 2;
t172 = t147 * t130;
t175 = t136 * t131;
t41 = -pkin(4) * t172 + pkin(7) * t175 + t47;
t14 = t138 * t41 - t141 * t185;
t15 = t185 * t138 + t141 * t41;
t8 = t138 * t15 - t141 * t14;
t183 = t134 * t8;
t59 = -t135 * t82 + t137 * t150;
t55 = -t131 * pkin(3) - t130 * qJ(4) + qJDD(4) - t59;
t43 = -pkin(4) * t175 + t55 + (-t130 * t145 - t172) * pkin(7);
t182 = t138 * t43;
t70 = qJDD(5) - t76;
t181 = t138 * t70;
t180 = t141 * t43;
t179 = t141 * t70;
t178 = t134 * t131;
t176 = t135 * t131;
t173 = t137 * t131;
t171 = t97 * qJD(5);
t170 = t99 * qJD(5);
t167 = qJD(5) * t138;
t166 = qJD(5) * t141;
t11 = t135 * t22 - t137 * t55;
t165 = pkin(2) * t11 - pkin(3) * t55 + qJ(4) * t22;
t109 = -t137 * t130 - t176;
t163 = pkin(2) * t109 - t60;
t9 = t138 * t14 + t141 * t15;
t103 = t186 * t187;
t86 = -t135 * t103 + t136 * t173;
t161 = pkin(2) * t86 + pkin(3) * t175 - qJ(4) * t103 - t136 * t55;
t63 = t188 * t131;
t94 = t153 * t131;
t56 = t138 * t63 - t141 * t94;
t57 = t138 * t94 + t141 * t63;
t31 = -t134 * t56 + t136 * t57;
t65 = -t95 - t96;
t17 = t135 * t31 - t137 * t65;
t160 = t134 * (-pkin(7) * t56 - t8) + t136 * (-pkin(4) * t65 + pkin(7) * t57 + t9) - pkin(3) * t65 + qJ(4) * t31 + pkin(2) * t17;
t144 = qJD(5) ^ 2;
t68 = -t144 - t95;
t53 = t138 * t68 + t189;
t54 = t141 * t68 - t190;
t29 = -t134 * t53 + t136 * t54;
t72 = -t63 + 0.2e1 * t170;
t19 = t135 * t29 - t137 * t72;
t159 = t134 * (-pkin(7) * t53 + t182) + t136 * (-pkin(4) * t72 + pkin(7) * t54 - t180) - pkin(3) * t72 + qJ(4) * t29 + pkin(2) * t19;
t91 = -t96 - t144;
t61 = t141 * t91 - t181;
t62 = -t138 * t91 - t179;
t37 = -t134 * t61 + t136 * t62;
t74 = t94 + 0.2e1 * t171;
t25 = t135 * t37 - t137 * t74;
t158 = t134 * (-pkin(7) * t61 + t180) + t136 * (-pkin(4) * t74 + pkin(7) * t62 + t182) - pkin(3) * t74 + qJ(4) * t37 + pkin(2) * t25;
t154 = t135 * t130 - t173;
t157 = -pkin(2) * t154 + t59;
t123 = t145 * t131;
t124 = t147 * t131;
t108 = t124 + t123;
t81 = t135 * t108 + t137 * t111;
t156 = pkin(2) * t81 + pkin(3) * t111 + qJ(4) * t108 + t22;
t102 = t134 * t111;
t85 = t135 * t102 - t134 * t173;
t152 = pkin(2) * t85 - pkin(3) * t178 + qJ(4) * t102 + t134 * t55;
t4 = t136 * t9 - t183;
t2 = t135 * t4 - t137 * t43;
t151 = pkin(2) * t2 - pkin(7) * t183 + qJ(4) * t4 - pkin(3) * t43 + t136 * (-pkin(4) * t43 + pkin(7) * t9);
t114 = 0.2e1 * t134 * t175;
t90 = -t96 + t144;
t89 = t95 - t144;
t75 = t94 + t171;
t73 = t63 - t170;
t48 = (t134 * (t138 * t99 + t141 * t97) + t136 * (t138 * t97 - t141 * t99)) * qJD(5);
t39 = t134 * (t141 * t75 - t99 * t167) + t136 * (t138 * t75 + t99 * t166);
t38 = t134 * (-t138 * t73 - t97 * t166) + t136 * (t141 * t73 - t97 * t167);
t36 = t134 * (-t138 * t90 + t189) + t136 * (t141 * t90 + t190);
t35 = t134 * (t141 * t89 - t181) + t136 * (t138 * t89 + t179);
t33 = t135 * t60 + t137 * t59;
t32 = pkin(2) * t33;
t30 = t134 * (-t138 * t74 - t141 * t72) + t136 * (-t138 * t72 + t141 * t74);
t1 = [0, 0, 0, 0, 0, qJDD(1), t162, t155, 0, 0, 0, 0, 0, 0, 0, t131, pkin(1) * (-t139 * t130 + t142 * t131) + t87, pkin(1) * (-t142 * t130 - t139 * t131) - t88, 0, pkin(1) * (t139 * t88 + t142 * t87), 0, 0, 0, 0, 0, t131, pkin(1) * (t139 * t109 - t142 * t154) + t157, pkin(1) * (t142 * t109 + t139 * t154) + t163, 0, pkin(1) * (t139 * (-t135 * t59 + t137 * t60) + t142 * t33) + t32, t123, t114, 0, t124, 0, 0, pkin(1) * (t139 * (-t137 * t103 - t135 * t175) + t142 * t86) + t161, pkin(1) * (t139 * (t137 * t102 + t134 * t176) + t142 * t85) + t152, pkin(1) * (t139 * (t137 * t108 - t135 * t111) + t142 * t81) + t156, pkin(1) * (t139 * (t135 * t55 + t137 * t22) + t142 * t11) + t165, t39, t30, t36, t38, t35, t48, pkin(1) * (t139 * (t135 * t72 + t137 * t29) + t142 * t19) + t159, pkin(1) * (t139 * (t135 * t74 + t137 * t37) + t142 * t25) + t158, pkin(1) * (t139 * (t135 * t65 + t137 * t31) + t142 * t17) + t160, pkin(1) * (t139 * (t135 * t43 + t137 * t4) + t142 * t2) + t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t131, t87, -t88, 0, 0, 0, 0, 0, 0, 0, t131, t157, t163, 0, t32, t123, t114, 0, t124, 0, 0, t161, t152, t156, t165, t39, t30, t36, t38, t35, t48, t159, t158, t160, t151; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t133, 0, 0, 0, 0, 0, 0, 0, 0, 0, t134 * t47 - t136 * t46, 0, 0, 0, 0, 0, 0, t134 * t54 + t136 * t53, t134 * t62 + t136 * t61, t134 * t57 + t136 * t56, t134 * t9 + t136 * t8; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t175, t178, -t111, t55, 0, 0, 0, 0, 0, 0, t72, t74, t65, t43; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t76, t96 - t95, t94, t76, t63, qJDD(5), -t14, -t15, 0, 0;];
tauJ_reg = t1;
