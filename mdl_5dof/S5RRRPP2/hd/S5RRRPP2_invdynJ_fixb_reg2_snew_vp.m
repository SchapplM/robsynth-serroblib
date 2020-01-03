% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RRRPP2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RRRPP2_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP2_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:51:59
% EndTime: 2019-12-31 20:52:02
% DurationCPUTime: 0.94s
% Computational Cost: add. (2499->185), mult. (3443->186), div. (0->0), fcn. (1609->6), ass. (0->115)
t113 = sin(qJ(2));
t116 = cos(qJ(2));
t112 = sin(qJ(3));
t115 = cos(qJ(3));
t105 = qJD(1) + qJD(2);
t103 = t105 ^ 2;
t149 = t115 * t103 * t112;
t84 = qJDD(3) - t149;
t168 = t115 * t84;
t118 = qJD(3) ^ 2;
t108 = t112 ^ 2;
t162 = t108 * t103;
t88 = t118 + t162;
t54 = -t112 * t88 + t168;
t155 = qJD(3) * t115;
t146 = t105 * t155;
t104 = qJDD(1) + qJDD(2);
t159 = t112 * t104;
t68 = 0.2e1 * t146 + t159;
t28 = pkin(1) * (t113 * t54 + t116 * t68);
t176 = pkin(2) * t68 + pkin(7) * t54;
t194 = 2 * qJD(4);
t83 = qJDD(3) + t149;
t193 = pkin(4) * t83;
t160 = t112 * qJ(4);
t133 = -t115 * pkin(3) - t160;
t67 = t133 * t105;
t173 = t105 * t67;
t114 = sin(qJ(1));
t117 = cos(qJ(1));
t135 = t114 * g(1) - t117 * g(2);
t125 = qJDD(1) * pkin(1) + t135;
t136 = t117 * g(1) + t114 * g(2);
t82 = -qJD(1) ^ 2 * pkin(1) - t136;
t43 = t113 * t125 + t116 * t82;
t37 = -t103 * pkin(2) + t104 * pkin(7) + t43;
t31 = -t112 * g(3) + t115 * t37;
t134 = qJDD(3) * qJ(4) + qJD(3) * t194 + t115 * t173 + t31;
t109 = t115 ^ 2;
t157 = t108 + t109;
t77 = t157 * t103;
t192 = t134 - (t118 - t77) * pkin(3);
t74 = t157 * t104;
t40 = pkin(1) * (t113 * t74 + t116 * t77);
t179 = t115 * g(3);
t30 = t112 * t37 + t179;
t18 = t112 * t30 + t115 * t31;
t163 = t105 * t112;
t128 = -qJD(3) * pkin(4) - qJ(5) * t163;
t158 = t115 * t104;
t156 = qJD(3) * t105;
t98 = t112 * t156;
t69 = -t98 + t158;
t190 = t69 * qJ(5) - qJD(3) * t128;
t189 = pkin(3) * t88 + qJ(4) * t84;
t161 = t109 * t103;
t91 = -t118 - t161;
t188 = pkin(3) * t83 + qJ(4) * t91;
t49 = t168 + t112 * (-t118 + t161);
t175 = pkin(2) * t77 + pkin(7) * t74;
t187 = t128 * t163 + qJDD(5);
t186 = t69 * pkin(4) + t187;
t141 = t37 + t173;
t145 = qJDD(3) * pkin(3) + t118 * qJ(4) - qJDD(4);
t24 = t141 * t112 - t145 + t179;
t185 = pkin(3) + pkin(4);
t126 = t146 + t159;
t42 = -t113 * t82 + t116 * t125;
t36 = -t104 * pkin(2) - t103 * pkin(7) - t42;
t124 = -t69 * pkin(3) + t36 + (-t126 - t146) * qJ(4);
t20 = (pkin(3) * qJD(3) - (2 * qJD(4))) * t163 + t124;
t23 = -t118 * pkin(3) + t134;
t9 = t112 * t24 + t115 * t23;
t180 = -pkin(2) * t20 + pkin(7) * t9;
t178 = -pkin(2) * t36 + pkin(7) * t18;
t172 = t112 * t83;
t53 = t115 * t91 - t172;
t70 = -0.2e1 * t98 + t158;
t177 = pkin(2) * t70 + pkin(7) * t53;
t174 = qJ(4) * t77;
t169 = t115 * t68;
t165 = qJ(4) * t115;
t164 = qJ(5) * t104;
t154 = qJD(5) * t105;
t152 = t112 * t36 - t176;
t151 = -t115 * t36 + t177;
t150 = 0.2e1 * t154;
t12 = qJ(5) * t161 - t186 + t20;
t122 = -0.2e1 * t115 * t154 - t190 + t23;
t13 = -pkin(4) * t161 + t122;
t129 = -t126 * qJ(5) - t145 - t193;
t121 = t129 + t30;
t147 = qJ(5) * t155;
t14 = (t147 + (-0.2e1 * qJD(5) + t67) * t112) * t105 + t121;
t4 = t112 * t14 + t115 * t13;
t148 = t112 * (-qJ(4) * t12 - qJ(5) * t14) + t115 * (-qJ(5) * t13 - t185 * t12) - pkin(2) * t12 + pkin(7) * t4;
t95 = t112 * t150;
t144 = t112 * (-t174 + t95 + (-qJ(5) * t156 - g(3)) * t115 + (-t141 + t164) * t112 - t129) + t115 * ((t150 + t164) * t115 + (-t77 + t161) * pkin(4) + t190 - t192) - t175;
t123 = t163 * t194 - t124;
t120 = (t70 - t98) * pkin(3) + t123;
t143 = t112 * (qJ(4) * t70 + qJ(5) * t83) + t115 * ((-t91 - t161) * qJ(5) + (t69 + t70) * pkin(4) + t120 + t187) + t177;
t119 = -pkin(3) * t98 + qJ(4) * t68 + t123;
t138 = t88 - t161;
t142 = t112 * (t138 * qJ(5) + t119 + t186) + t115 * (-qJ(5) * t84 + t185 * t68) + t176;
t140 = t112 * (t24 + t174) + t115 * t192 + t175;
t139 = t175 + t18;
t38 = t112 * t70 + t169;
t51 = -t115 * (-t118 + t162) + t172;
t130 = pkin(3) * t169 + t112 * t119 + t176;
t127 = t115 * t120 + t70 * t160 + t177;
t78 = (t108 - t109) * t103;
t45 = t68 * t112;
t44 = (t69 - t98) * t115;
t27 = pkin(1) * (t113 * t53 + t116 * t70);
t1 = [0, 0, 0, 0, 0, qJDD(1), t135, t136, 0, 0, 0, 0, 0, 0, 0, t104, pkin(1) * (-t113 * t103 + t116 * t104) + t42, pkin(1) * (-t116 * t103 - t113 * t104) - t43, 0, pkin(1) * (t113 * t43 + t116 * t42), t45, t38, t51, t44, t49, 0, t27 + t151, t152 - t28, t40 + t139, pkin(1) * (t113 * t18 - t116 * t36) + t178, t45, t51, -t38, 0, -t49, t44, t127 + t27, t40 + t140, t130 + t28, pkin(1) * t113 * t9 + (-pkin(1) * t116 + t133) * t20 + t180, t45, -t38, -t51, t44, t49, 0, t27 + t143, t28 + t142, t144 - t40, pkin(1) * (t113 * t4 - t116 * t12) + t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t104, t42, -t43, 0, 0, t45, t38, t51, t44, t49, 0, t151, t152, t139, t178, t45, t51, -t38, 0, -t49, t44, t127, t140, t130, t133 * t20 + t180, t45, -t38, -t51, t44, t49, 0, t143, t142, t144, t148; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t149, t78, t159, t149, t158, qJDD(3), -t30, -t31, 0, 0, -t149, t159, -t78, qJDD(3), -t158, t149, -t24 + t188, (-pkin(3) * t112 + t165) * t104, t23 + t189, -pkin(3) * t24 + qJ(4) * t23, -t149, -t78, -t159, t149, t158, qJDD(3), t193 + t95 + (-t112 * t67 - t147) * t105 - t121 + t188, t138 * pkin(4) + t122 + t189, (t185 * t112 - t165) * t104, qJ(4) * t13 - t185 * t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t83, t159, -t88, t24, 0, 0, 0, 0, 0, 0, -t83, -t88, -t159, t14; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t70, t68, -t77, -t12;];
tauJ_reg = t1;
