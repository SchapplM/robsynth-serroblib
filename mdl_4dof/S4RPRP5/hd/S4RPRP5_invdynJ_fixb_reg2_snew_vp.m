% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% qJDD [4x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:45:01
% EndTime: 2019-12-31 16:45:05
% DurationCPUTime: 1.11s
% Computational Cost: add. (1466->174), mult. (3692->212), div. (0->0), fcn. (2431->6), ass. (0->108)
t94 = sin(pkin(6));
t96 = sin(qJ(3));
t138 = t94 * t96;
t95 = cos(pkin(6));
t98 = cos(qJ(3));
t77 = (-t95 * t98 + t138) * qJD(1);
t137 = t94 * t98;
t109 = t95 * t96 + t137;
t79 = t109 * qJD(1);
t139 = t79 * t77;
t148 = qJDD(3) + t139;
t135 = t96 * t148;
t74 = t79 ^ 2;
t99 = qJD(3) ^ 2;
t152 = -t74 - t99;
t17 = -t152 * t98 + t135;
t129 = t98 * t148;
t19 = t152 * t96 + t129;
t177 = qJ(2) * (t94 * t17 - t95 * t19);
t176 = pkin(5) * t17;
t175 = pkin(5) * t19;
t147 = t77 ^ 2;
t61 = t147 - t99;
t173 = t94 * (-t98 * t61 + t135) - t95 * (t96 * t61 + t129);
t149 = qJDD(3) - t139;
t134 = t96 * t149;
t151 = -t147 - t99;
t155 = t151 * t98 - t134;
t170 = pkin(5) * t155;
t37 = t98 * t149;
t156 = t151 * t96 + t37;
t169 = pkin(5) * t156;
t121 = t95 * qJDD(1);
t122 = t94 * qJDD(1);
t75 = -t98 * t121 + t96 * t122;
t76 = t109 * qJDD(1);
t158 = -t96 * t75 - t98 * t76;
t168 = pkin(5) * t158;
t157 = -t98 * t75 + t96 * t76;
t30 = t74 + t147;
t165 = pkin(2) * t30 + pkin(5) * t157;
t62 = -t74 + t99;
t164 = t94 * (-t96 * t62 + t37) + t95 * (t98 * t62 + t134);
t123 = t79 * qJD(3);
t49 = t75 + 0.2e1 * t123;
t163 = qJ(2) * (t155 * t95 - t156 * t94) - pkin(1) * t49;
t162 = pkin(1) * t30 + qJ(2) * (t157 * t95 - t158 * t94);
t100 = qJD(1) ^ 2;
t141 = cos(qJ(1));
t97 = sin(qJ(1));
t108 = t141 * g(1) + t97 * g(2);
t114 = -t100 * pkin(1) + qJDD(1) * qJ(2) + 0.2e1 * qJD(1) * qJD(2) - t108;
t70 = qJD(3) * t77;
t52 = t76 - t70;
t153 = -t70 + t52;
t159 = t153 * qJ(4);
t124 = t100 * qJ(2);
t125 = qJDD(1) * pkin(1);
t89 = t94 ^ 2;
t90 = t95 ^ 2;
t127 = t89 + t90;
t119 = t97 * g(1) - t141 * g(2);
t110 = -qJDD(2) + t119;
t71 = t110 + t124 + t125;
t154 = t127 * t124 - t125 - t71;
t150 = t74 - t147;
t146 = 2 * qJD(4);
t145 = pkin(3) * t98;
t50 = -t75 - t123;
t144 = t50 * pkin(3);
t140 = pkin(2) * t100;
t142 = t95 * g(3);
t107 = -t142 + (-pkin(5) * qJDD(1) + t95 * t140 - t114) * t94;
t113 = -t94 * g(3) + t114 * t95;
t34 = pkin(5) * t121 - t90 * t140 + t113;
t14 = -t98 * t107 + t96 * t34;
t15 = t96 * t107 + t98 * t34;
t7 = -t98 * t14 + t96 * t15;
t143 = t94 * t7;
t43 = (t95 * pkin(2) + pkin(1)) * qJDD(1) + (t127 * pkin(5) + qJ(2)) * t100 + t110;
t136 = t96 * t43;
t133 = t96 * t49;
t130 = t98 * t43;
t128 = t98 * t49;
t126 = qJ(4) * t98;
t118 = -qJ(4) * t96 - pkin(2);
t8 = t96 * t14 + t98 * t15;
t117 = t95 * t113 + t94 * (t114 * t94 + t142);
t39 = t77 * pkin(3) - t79 * qJ(4);
t112 = qJDD(3) * qJ(4) + qJD(3) * t146 - t77 * t39 + t15;
t11 = -qJDD(3) * pkin(3) - t99 * qJ(4) + t79 * t39 + qJDD(4) + t14;
t106 = t94 * (-t96 * t50 + t98 * t70) + t95 * (t98 * t50 + t96 * t70);
t59 = t96 * t123;
t105 = t94 * t59 + (-t77 * t137 + t95 * (-t77 * t96 - t79 * t98)) * qJD(3);
t104 = -pkin(3) * t123 + t79 * t146 + t43;
t103 = t104 + t159;
t86 = t90 * qJDD(1);
t85 = t89 * qJDD(1);
t81 = t127 * t100;
t51 = t76 - 0.2e1 * t70;
t12 = t94 * (t98 * t52 - t59) + t95 * (t98 * t123 + t96 * t52);
t10 = -t99 * pkin(3) + t112;
t9 = t103 + t144;
t6 = (-t49 + t50) * pkin(3) + t103;
t5 = qJ(4) * t30 + t11;
t4 = (t30 - t99) * pkin(3) + t112;
t3 = t104 + t144 + 0.2e1 * t159;
t1 = [0, 0, 0, 0, 0, qJDD(1), t119, t108, 0, 0, t85, 0.2e1 * t94 * t121, 0, t86, 0, 0, -t154 * t95, t154 * t94, pkin(1) * t81 + qJ(2) * (t86 + t85) + t117, pkin(1) * t71 + qJ(2) * t117, t12, t94 * (-t96 * t51 - t128) + t95 * (t98 * t51 - t133), t164, t106, -t173, t105, t94 * (-t136 - t169) + t95 * (-pkin(2) * t49 + t130 + t170) + t163, t94 * (-t130 + t176) + t95 * (-pkin(2) * t51 - t136 - t175) - pkin(1) * t51 + t177, t94 * (-t7 - t168) + t95 * (t165 + t8) + t162, -pkin(5) * t143 + t95 * (pkin(2) * t43 + pkin(5) * t8) + pkin(1) * t43 + qJ(2) * (t95 * t8 - t143), t12, t164, t94 * (t153 * t96 + t128) + t95 * (-t153 * t98 + t133), t105, t173, t106, t94 * (-t49 * t126 - t96 * t6 - t169) + t95 * (t118 * t49 + t98 * t6 + t170) + t163, t94 * (-t96 * t4 + t98 * t5 - t168) + t95 * (t98 * t4 + t96 * t5 + t165) + t162, t94 * (t98 * t3 - t176) + t95 * (t96 * t3 + t175) - t177 + (-pkin(3) * t138 + t95 * (pkin(2) + t145) + pkin(1)) * t153, (t94 * (-pkin(3) * t96 + t126) + t95 * (-t118 + t145) + pkin(1)) * t9 + (qJ(2) + pkin(5)) * (-t94 * (t96 * t10 - t98 * t11) + t95 * (t98 * t10 + t96 * t11)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t121, t122, -t81, -t71, 0, 0, 0, 0, 0, 0, t49, t51, -t30, -t43, 0, 0, 0, 0, 0, 0, t49, -t30, -t153, -t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t139, t150, t76, -t139, -t75, qJDD(3), -t14, -t15, 0, 0, t139, t70 + t52, -t150, qJDD(3), t75, -t139, pkin(3) * t149 + qJ(4) * t151 - t11, -pkin(3) * t76 - qJ(4) * t75, qJ(4) * t148 + (-t152 - t99) * pkin(3) + t112, -pkin(3) * t11 + qJ(4) * t10; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t76, t152, t11;];
tauJ_reg = t1;
