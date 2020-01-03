% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S4RRPP3
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
%   pkin=[a2,a3,a4,d1,d2,theta3]';
% 
% Output:
% tauJ_reg [4x(4*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:58
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S4RRPP3_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(6,1)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_snew_vp: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_snew_vp: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_snew_vp: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP3_invdynJ_fixb_reg2_snew_vp: pkin has to be [6x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:57:51
% EndTime: 2019-12-31 16:57:55
% DurationCPUTime: 1.51s
% Computational Cost: add. (1970->203), mult. (4705->241), div. (0->0), fcn. (2973->6), ass. (0->123)
t107 = sin(qJ(2));
t109 = cos(qJ(2));
t106 = cos(pkin(6));
t105 = sin(pkin(6));
t134 = qJD(1) * t107;
t79 = -t106 * t109 * qJD(1) + t105 * t134;
t136 = t107 * t106;
t81 = (t109 * t105 + t136) * qJD(1);
t155 = t81 * t79;
t163 = qJDD(2) + t155;
t150 = t105 * t163;
t110 = qJD(2) ^ 2;
t78 = t81 ^ 2;
t167 = -t78 - t110;
t21 = -t106 * t167 + t150;
t145 = t106 * t163;
t23 = t105 * t167 + t145;
t201 = pkin(5) * (t107 * t21 - t109 * t23);
t200 = pkin(2) * t21;
t199 = qJ(3) * t21;
t198 = qJ(3) * t23;
t162 = t79 ^ 2;
t64 = t162 - t110;
t196 = t107 * (-t106 * t64 + t150) - t109 * (t105 * t64 + t145);
t131 = qJD(1) * qJD(2);
t128 = t109 * t131;
t94 = t107 * qJDD(1);
t87 = t94 + t128;
t129 = t107 * t131;
t95 = t109 * qJDD(1);
t88 = t95 - t129;
t125 = t105 * t87 - t106 * t88;
t75 = qJD(2) * t81;
t39 = t125 - t75;
t140 = qJD(2) * t79;
t58 = t105 * t88 + t106 * t87;
t43 = t58 + t140;
t182 = -t105 * t39 - t106 * t43;
t195 = pkin(2) * t182;
t194 = qJ(3) * t182;
t181 = t105 * t43 - t106 * t39;
t37 = t78 + t162;
t191 = pkin(2) * t37 + qJ(3) * t181;
t190 = pkin(1) * t37 + pkin(5) * (-t107 * t182 + t109 * t181);
t164 = qJDD(2) - t155;
t149 = t105 * t164;
t165 = -t162 - t110;
t170 = t106 * t165 - t149;
t185 = qJ(3) * t170;
t44 = t106 * t164;
t172 = t105 * t165 + t44;
t184 = qJ(3) * t172;
t138 = qJD(3) * t81;
t183 = pkin(2) * t172 - 0.2e1 * t138;
t168 = t58 - t140;
t38 = t125 + t75;
t180 = t107 * (t105 * t168 + t106 * t38) - t109 * (-t105 * t38 + t106 * t168);
t65 = -t78 + t110;
t179 = t107 * (-t105 * t65 + t44) + t109 * (t106 * t65 + t149);
t178 = pkin(5) * (-t107 * t172 + t109 * t170) - pkin(1) * t38;
t173 = t168 * qJ(4);
t166 = t78 - t162;
t111 = qJD(1) ^ 2;
t135 = t107 * t111;
t108 = sin(qJ(1));
t159 = cos(qJ(1));
t120 = t159 * g(1) + t108 * g(2);
t137 = qJDD(1) * pkin(5);
t84 = -t111 * pkin(1) - t120 + t137;
t144 = t107 * t84;
t114 = qJDD(2) * pkin(2) - t87 * qJ(3) - t144 + (pkin(2) * t135 + qJ(3) * t131 - g(3)) * t109;
t118 = qJD(2) * pkin(2) - qJ(3) * t134;
t61 = -t107 * g(3) + t109 * t84;
t103 = t109 ^ 2;
t97 = t103 * t111;
t31 = -pkin(2) * t97 + t88 * qJ(3) - qJD(2) * t118 + t61;
t14 = -0.2e1 * qJD(3) * t79 + t105 * t114 + t106 * t31;
t161 = 2 * qJD(4);
t160 = t125 * pkin(3);
t158 = pkin(3) * t105;
t157 = pkin(3) * t106;
t126 = t105 * t31 - t106 * t114;
t13 = t126 + 0.2e1 * t138;
t3 = t105 * t14 - t106 * t13;
t156 = t107 * t3;
t130 = t108 * g(1) - t159 * g(2);
t119 = qJDD(1) * pkin(1) + t130;
t35 = t88 * pkin(2) - qJDD(3) - t118 * t134 + (qJ(3) * t103 + pkin(5)) * t111 + t119;
t153 = t105 * t35;
t147 = t106 * t35;
t93 = t109 * t135;
t143 = t107 * (qJDD(2) + t93);
t142 = t109 * (qJDD(2) - t93);
t141 = qJ(4) * t106;
t133 = qJD(2) * t105;
t132 = qJD(2) * t106;
t127 = -qJ(4) * t105 - pkin(2);
t4 = t105 * t13 + t106 * t14;
t60 = t109 * g(3) + t144;
t124 = t107 * t60 + t109 * t61;
t45 = t79 * pkin(3) - t81 * qJ(4);
t121 = qJDD(2) * qJ(4) + qJD(2) * t161 - t79 * t45 + t14;
t117 = -qJDD(2) * pkin(3) - t110 * qJ(4) + qJDD(4) + t126;
t10 = (0.2e1 * qJD(3) + t45) * t81 + t117;
t116 = t107 * (t105 * t125 + t79 * t132) + t109 * (-t106 * t125 + t79 * t133);
t63 = t81 * t133;
t115 = t107 * t63 + (-t79 * t136 + t109 * (-t105 * t79 - t106 * t81)) * qJD(2);
t113 = -pkin(3) * t75 + t81 * t161 + t35;
t112 = t113 + t173;
t102 = t107 ^ 2;
t96 = t102 * t111;
t89 = t95 - 0.2e1 * t129;
t86 = t94 + 0.2e1 * t128;
t83 = t111 * pkin(5) + t119;
t15 = t107 * (t106 * t58 - t63) + t109 * (t105 * t58 + t81 * t132);
t11 = t112 - t160;
t9 = -t110 * pkin(3) + t121;
t8 = (-t38 - t125) * pkin(3) + t112;
t7 = t113 - t160 + 0.2e1 * t173;
t6 = qJ(4) * t37 + t10;
t5 = (-t110 + t37) * pkin(3) + t121;
t1 = -t106 * t10 + t105 * t9;
t2 = [0, 0, 0, 0, 0, qJDD(1), t130, t120, 0, 0, (t87 + t128) * t107, t107 * t89 + t109 * t86, t143 + t109 * (-t96 + t110), (t88 - t129) * t109, t107 * (t97 - t110) + t142, 0, t109 * t83 + pkin(1) * t89 + pkin(5) * (t109 * (-t97 - t110) - t143), -t107 * t83 - pkin(1) * t86 + pkin(5) * (-t142 - t107 * (-t96 - t110)), pkin(1) * (t96 + t97) + (t102 + t103) * t137 + t124, pkin(1) * t83 + pkin(5) * t124, t15, -t180, t179, t116, -t196, t115, t107 * (-t153 - t184) + t109 * (-pkin(2) * t38 + t147 + t185) + t178, t107 * (-t147 + t199) + t109 * (-pkin(2) * t168 - t153 - t198) - pkin(1) * t168 + t201, t107 * (-t3 - t194) + t109 * (t191 + t4) + t190, -qJ(3) * t156 + t109 * (pkin(2) * t35 + qJ(3) * t4) + pkin(1) * t35 + pkin(5) * (t109 * t4 - t156), t15, t179, t180, t115, t196, t116, t107 * (-t105 * t8 - t38 * t141 - t184) + t109 * (t106 * t8 + t127 * t38 + t185) + t178, t107 * (-t105 * t5 + t106 * t6 - t194) + t109 * (t105 * t6 + t106 * t5 + t191) + t190, t107 * (t106 * t7 - t199) + t109 * (t105 * t7 + t198) - t201 + (-t107 * t158 + t109 * (pkin(2) + t157) + pkin(1)) * t168, (t107 * (t141 - t158) + t109 * (-t127 + t157) + pkin(1)) * t11 + (pkin(5) + qJ(3)) * (-t107 * t1 + t109 * (t105 * t10 + t106 * t9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t93, t96 - t97, t94, t93, t95, qJDD(2), -t60, -t61, 0, 0, t155, t166, t43, -t155, -t39, qJDD(2), -t126 + t183, -t14 - t200, t195, pkin(2) * t3, t155, t43, -t166, qJDD(2), t39, -t155, pkin(3) * t164 + qJ(4) * t165 - t81 * t45 - t117 + t183, -pkin(3) * t43 - qJ(4) * t39 + t195, t200 + qJ(4) * t163 + (-t110 - t167) * pkin(3) + t121, pkin(2) * t1 - pkin(3) * t10 + qJ(4) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t38, t168, -t37, -t35, 0, 0, 0, 0, 0, 0, t38, -t37, -t168, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t164, t43, t167, t10;];
tauJ_reg = t2;
