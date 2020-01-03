% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPRRP12
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:57
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPRRP12_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRRP12_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:57:20
% EndTime: 2019-12-31 18:57:25
% DurationCPUTime: 1.30s
% Computational Cost: add. (3233->206), mult. (6335->247), div. (0->0), fcn. (3813->6), ass. (0->156)
t121 = sin(qJ(4));
t124 = cos(qJ(4));
t125 = cos(qJ(3));
t160 = qJD(1) * t125;
t101 = -t124 * qJD(3) + t121 * t160;
t114 = t125 * qJDD(1);
t122 = sin(qJ(3));
t156 = qJD(1) * qJD(3);
t151 = t122 * t156;
t106 = t114 - t151;
t73 = -t101 * qJD(4) + t121 * qJDD(3) + t124 * t106;
t112 = t122 * qJD(1) + qJD(4);
t91 = t112 * t101;
t59 = t73 + t91;
t197 = qJ(5) * t59;
t188 = pkin(6) + pkin(1);
t150 = t125 * t156;
t154 = t122 * qJDD(1);
t105 = -t150 - t154;
t100 = qJDD(4) - t105;
t103 = t121 * qJD(3) + t124 * t160;
t80 = t103 * t101;
t190 = t100 - t80;
t196 = pkin(4) * t190;
t195 = t121 * t190;
t194 = t124 * t190;
t155 = qJD(2) * qJD(1);
t116 = 0.2e1 * t155;
t118 = qJDD(1) * qJ(2);
t123 = sin(qJ(1));
t126 = cos(qJ(1));
t143 = t126 * g(1) + t123 * g(2);
t137 = -t118 + t143;
t134 = t116 - t137;
t140 = -t106 + t151;
t141 = -t105 + t150;
t128 = qJD(1) ^ 2;
t189 = t188 * t128;
t49 = t141 * pkin(3) + t140 * pkin(7) + t134 - t189;
t127 = qJD(3) ^ 2;
t180 = pkin(3) * t122;
t144 = -pkin(7) * t125 + t180;
t135 = t128 * t144;
t149 = t123 * g(1) - t126 * g(2);
t142 = qJDD(2) - t149;
t162 = t128 * qJ(2);
t132 = t142 - t162;
t89 = -t188 * qJDD(1) + t132;
t76 = t125 * g(3) - t122 * t89;
t62 = -t127 * pkin(3) + qJDD(3) * pkin(7) - t122 * t135 - t76;
t26 = t121 * t49 + t124 * t62;
t148 = -t124 * qJDD(3) + t121 * t106;
t72 = -t103 * qJD(4) - t148;
t83 = t112 * pkin(4) - t103 * qJ(5);
t136 = t72 * qJ(5) - 0.2e1 * qJD(5) * t101 - t112 * t83 + t26;
t111 = t112 ^ 2;
t99 = t103 ^ 2;
t77 = -t99 - t111;
t98 = t101 ^ 2;
t176 = t77 + t98;
t193 = t176 * pkin(4) - t136;
t191 = t73 - t91;
t55 = (qJD(4) - t112) * t103 + t148;
t74 = -t111 - t98;
t37 = t121 * t74 + t194;
t187 = pkin(3) * t37;
t67 = t100 + t80;
t174 = t121 * t67;
t41 = t124 * t77 - t174;
t186 = pkin(3) * t41;
t25 = t121 * t62 - t124 * t49;
t133 = t196 - t25 - t197;
t158 = qJD(5) * t103;
t95 = -0.2e1 * t158;
t10 = t133 + t95;
t185 = pkin(4) * t10;
t184 = pkin(4) * t59;
t30 = -t121 * t55 - t124 * t59;
t183 = pkin(7) * t30;
t182 = pkin(7) * t37;
t181 = pkin(7) * t41;
t31 = t121 * t59 - t124 * t55;
t65 = -t98 - t99;
t179 = -pkin(3) * t65 + pkin(7) * t31;
t38 = t124 * t74 - t195;
t54 = (qJD(4) + t112) * t103 + t148;
t178 = -pkin(3) * t54 + pkin(7) * t38;
t172 = t124 * t67;
t42 = -t121 * t77 - t172;
t177 = -pkin(3) * t191 + pkin(7) * t42;
t75 = t122 * g(3) + t125 * t89;
t61 = qJDD(3) * pkin(3) + t127 * pkin(7) - t125 * t135 + t75;
t175 = t121 * t61;
t173 = t124 * t61;
t171 = qJ(5) * t121;
t170 = qJ(5) * t124;
t169 = qJDD(1) * pkin(1);
t168 = t112 * t121;
t167 = t112 * t124;
t119 = t122 ^ 2;
t166 = t119 * t128;
t120 = t125 ^ 2;
t165 = t120 * t128;
t152 = t122 * t128 * t125;
t164 = t122 * (qJDD(3) + t152);
t163 = t125 * (qJDD(3) - t152);
t161 = t119 + t120;
t153 = t122 * t80;
t9 = t121 * t25 + t124 * t26;
t15 = t122 * t31 - t125 * t65;
t147 = qJ(2) * t30 - t188 * t15;
t18 = t122 * t38 - t125 * t54;
t146 = qJ(2) * t37 - t188 * t18;
t21 = t122 * t42 - t125 * t191;
t145 = qJ(2) * t41 - t188 * t21;
t139 = t121 * t26 - t124 * t25;
t48 = -t122 * t76 + t125 * t75;
t131 = t133 + t196;
t130 = -t103 * t83 - qJDD(5) + t61;
t129 = -t72 * pkin(4) - t130;
t108 = t161 * qJDD(1);
t107 = t114 - 0.2e1 * t151;
t104 = 0.2e1 * t150 + t154;
t96 = 0.2e1 * t158;
t92 = -t132 + t169;
t87 = -t99 + t111;
t86 = t98 - t111;
t85 = t137 - 0.2e1 * t155 + t189;
t82 = -t164 + t125 * (-t127 - t165);
t81 = t122 * (-t127 - t166) + t163;
t78 = t99 - t98;
t63 = (-t101 * t121 - t103 * t124) * t112;
t51 = t103 * t167 + t121 * t73;
t50 = t101 * t168 + t124 * t72;
t45 = t122 * t100 + t125 * (-t101 * t124 + t103 * t121) * t112;
t44 = t121 * t86 + t172;
t43 = t124 * t87 + t195;
t34 = t125 * (-t103 * t168 + t124 * t73) + t153;
t33 = t125 * (t101 * t167 - t121 * t72) - t153;
t32 = -pkin(4) * t191 - qJ(5) * t67;
t29 = -t121 * t54 + t124 * t191;
t23 = t125 * (t124 * t86 - t174) - t122 * t55;
t22 = t125 * (-t121 * t87 + t194) + t122 * t59;
t19 = t98 * qJ(5) - t129;
t16 = t125 * (-t121 * t191 - t124 * t54) + t122 * t78;
t13 = -t176 * qJ(5) + t129;
t12 = -t98 * pkin(4) + t136;
t11 = (t74 + t98) * qJ(5) + (-t54 + t72) * pkin(4) + t130;
t7 = -t133 + t96 + t197;
t6 = t122 * t9 + t125 * t61;
t5 = -qJ(5) * t55 + (-t65 - t98) * pkin(4) + t136;
t4 = pkin(4) * t19 + qJ(5) * t12;
t3 = -t121 * t10 + t124 * t12;
t2 = t124 * t10 + t121 * t12;
t1 = t122 * t3 + t125 * t19;
t8 = [0, 0, 0, 0, 0, qJDD(1), t149, t143, 0, 0, qJDD(1), 0, 0, 0, 0, 0, 0, t142 - 0.2e1 * t169, t116 + 0.2e1 * t118 - t143, pkin(1) * t92 + qJ(2) * (-t128 * pkin(1) + t134), -t140 * t125, -t125 * t104 - t122 * t107, t163 - t122 * (t127 - t165), t141 * t122, t125 * (-t127 + t166) - t164, 0, qJ(2) * t104 - t122 * t85 - t188 * t81, qJ(2) * t107 - t125 * t85 - t188 * t82, t188 * t108 - t161 * t162 - t48, -qJ(2) * t85 - t188 * t48, t34, t16, t22, t33, t23, t45, t125 * (-t175 - t182) - t122 * (t25 - t187) + t146, t125 * (-t173 - t181) - t122 * (t26 - t186) + t145, t125 * (-t139 - t183) + t30 * t180 + t147, -t188 * t6 + (qJ(2) + t144) * t139, t34, t16, t22, t33, t23, t45, t125 * (-t121 * t11 - t170 * t190 - t182) - t122 * (-t131 + t96 - t187) + t146, t125 * (-t121 * t32 + t124 * t13 - t181) - t122 * (-t186 - t193) + t145, t125 * (-t121 * t5 + t124 * t7 - t183) - t122 * (-pkin(3) * t30 + t184) + t147, t125 * (-pkin(7) * t2 - t10 * t170 - t121 * t4) - t122 * (-pkin(3) * t2 - t185) + qJ(2) * t2 - t188 * t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, qJDD(1), -t128, -t92, 0, 0, 0, 0, 0, 0, t81, t82, -t108, t48, 0, 0, 0, 0, 0, 0, t18, t21, t15, t6, 0, 0, 0, 0, 0, 0, t18, t21, t15, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t152, (-t119 + t120) * t128, t114, -t152, -t154, qJDD(3), t75, t76, 0, 0, t51, t29, t43, t50, t44, t63, t173 + t178, -t175 + t177, t9 + t179, pkin(3) * t61 + pkin(7) * t9, t51, t29, t43, t50, t44, t63, t124 * t11 - t171 * t190 + t178, t121 * t13 + t124 * t32 + t177, t121 * t7 + t124 * t5 + t179, pkin(3) * t19 + pkin(7) * t3 - t10 * t171 + t124 * t4; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t80, t78, t59, -t80, -t55, t100, -t25, -t26, 0, 0, t80, t78, t59, -t80, -t55, t100, t131 + t95, t193, -t184, t185; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t54, t191, t65, -t19;];
tauJ_reg = t8;
