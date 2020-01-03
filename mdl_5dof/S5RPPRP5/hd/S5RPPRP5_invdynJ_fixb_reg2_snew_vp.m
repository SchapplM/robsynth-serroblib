% Calculate inertial parameters regressor of inverse dynamics joint torque vector with Newton-Euler for
% S5RPPRP5
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
%   pkin=[a2,a3,a4,a5,d1,d4,theta2]';
% 
% Output:
% tauJ_reg [5x(5*10)]
%   inertial parameter regressor of inverse dynamics joint torque vector

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauJ_reg = S5RPPRP5_invdynJ_fixb_reg2_snew_vp(qJ, qJD, qJDD, g, pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_snew_vp: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_snew_vp: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_snew_vp: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_snew_vp: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRP5_invdynJ_fixb_reg2_snew_vp: pkin has to be [7x1] (double)');

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_tauJ_regressor_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:53:43
% EndTime: 2019-12-31 17:53:49
% DurationCPUTime: 1.92s
% Computational Cost: add. (2043->207), mult. (5088->240), div. (0->0), fcn. (3201->6), ass. (0->137)
t104 = sin(pkin(7));
t105 = cos(pkin(7));
t108 = cos(qJ(4));
t106 = sin(qJ(4));
t127 = t104 * t106 + t105 * t108;
t78 = t127 * qJD(1);
t150 = qJD(1) * t105;
t151 = qJD(1) * t104;
t80 = -t106 * t150 + t108 * t151;
t176 = t80 * t78;
t203 = qJDD(4) + t176;
t170 = t106 * t203;
t109 = qJD(4) ^ 2;
t75 = t80 ^ 2;
t191 = -t75 - t109;
t16 = -t108 * t191 + t170;
t164 = t108 * t203;
t18 = t106 * t191 + t164;
t220 = qJ(2) * (t104 * t16 + t105 * t18);
t219 = pkin(6) * t16;
t218 = pkin(6) * t18;
t185 = t78 ^ 2;
t60 = t185 - t109;
t216 = t104 * (-t108 * t60 + t170) + t105 * (t106 * t60 + t164);
t204 = qJDD(4) - t176;
t163 = t108 * t204;
t169 = t106 * t204;
t61 = -t75 + t109;
t213 = t104 * (-t106 * t61 + t163) + t105 * (-t108 * t61 - t169);
t159 = qJD(4) * t78;
t147 = t105 * qJDD(1);
t97 = t104 * qJDD(1);
t77 = -t106 * t147 + t108 * t97;
t121 = t77 - t159;
t32 = t121 - t159;
t207 = qJ(5) * t32;
t107 = sin(qJ(1));
t181 = cos(qJ(1));
t140 = t107 * g(1) - t181 * g(2);
t110 = qJD(1) ^ 2;
t155 = t110 * qJ(2);
t122 = -qJDD(2) + t140 + t155;
t119 = 0.2e1 * qJD(3) * t151 + t122;
t157 = t104 * qJ(3);
t182 = pkin(2) + pkin(3);
t123 = t105 * t182 + pkin(1) + t157;
t111 = t104 ^ 2;
t153 = t111 * t110;
t113 = t105 ^ 2;
t154 = t110 * t113;
t29 = t123 * qJDD(1) + t119 + (-t153 - t154) * pkin(6);
t152 = t80 * qJD(4);
t187 = t127 * qJDD(1);
t51 = -t187 - t152;
t212 = -t51 * pkin(4) - t207 + t29;
t188 = -t185 - t109;
t194 = t108 * t188 - t169;
t211 = pkin(6) * t194;
t195 = t106 * t77 - t108 * t187;
t210 = pkin(6) * t195;
t196 = t106 * t188 + t163;
t209 = pkin(6) * t196;
t197 = -t106 * t187 - t108 * t77;
t208 = pkin(6) * t197;
t50 = t187 + 0.2e1 * t152;
t202 = pkin(1) * t50 + qJ(2) * (t104 * t196 + t105 * t194);
t35 = t75 + t185;
t201 = qJ(2) * (t104 * t197 + t105 * t195) - t123 * t35;
t184 = 2 * qJD(2);
t125 = t181 * g(1) + t107 * g(2);
t82 = -t110 * pkin(1) + qJDD(1) * qJ(2) - t125;
t136 = qJD(1) * t184 + t82;
t193 = qJ(2) - pkin(6);
t190 = t75 - t185;
t189 = t111 + t113;
t87 = t189 * t110;
t183 = 2 * qJD(5);
t180 = pkin(4) * t106;
t179 = pkin(4) * t108;
t178 = t105 * pkin(2);
t177 = t105 * g(3);
t130 = -t157 - t178;
t85 = t130 * qJD(1);
t128 = t82 + (t184 + t85) * qJD(1);
t138 = qJDD(3) + t177;
t156 = t105 * t110;
t116 = (-pkin(3) * t156 - pkin(6) * qJDD(1) + t128) * t104 + t138;
t175 = t136 * t105;
t134 = -t104 * g(3) + t175;
t71 = t85 * t150;
t129 = t134 + t71;
t28 = -pkin(3) * t154 - pkin(6) * t147 + t129;
t13 = t106 * t116 + t108 * t28;
t96 = t111 * qJDD(1);
t98 = t113 * qJDD(1);
t174 = qJ(2) * (t98 + t96) + pkin(1) * t87;
t173 = -t189 * t105 * t155 + pkin(1) * t147;
t172 = qJ(2) * t104 * t87;
t171 = t106 * t29;
t168 = t106 * t50;
t165 = t108 * t29;
t162 = t108 * t50;
t158 = qJDD(1) * pkin(1);
t149 = qJD(4) * t108;
t148 = t106 * qJD(4);
t145 = t78 * t148;
t144 = t78 * t149;
t142 = t104 * t147;
t139 = t104 * (t136 * t104 + t177) + t105 * t134;
t12 = t106 * t28 - t108 * t116;
t137 = -qJ(5) * t108 + qJ(3);
t135 = qJ(5) * t106 + t182;
t39 = t78 * pkin(4) - t80 * qJ(5);
t133 = qJDD(4) * qJ(5) + qJD(4) * t183 - t78 * t39 + t13;
t3 = t106 * t13 - t108 * t12;
t4 = t106 * t12 + t108 * t13;
t126 = pkin(1) - t130;
t120 = t126 * qJDD(1);
t10 = -qJDD(4) * pkin(4) - t109 * qJ(5) + t80 * t39 + qJDD(5) + t12;
t58 = t80 * t148;
t59 = t80 * t149;
t118 = t104 * (t58 - t144) + t105 * (t59 + t145);
t117 = t104 * (-t106 * t51 + t144) + t105 * (-t108 * t51 - t145);
t115 = t80 * t183 - t212;
t72 = t122 + t158;
t52 = t77 - 0.2e1 * t159;
t40 = t120 + t119;
t33 = t128 * t104 + t138;
t14 = t104 * (t108 * t121 - t58) + t105 * (-t106 * t121 - t59);
t11 = (pkin(4) * qJD(4) - (2 * qJD(5))) * t80 + t212;
t9 = -t109 * pkin(4) + t133;
t8 = t115 + (-t50 - t152) * pkin(4);
t7 = -pkin(4) * t152 + t115 + t207;
t6 = qJ(5) * t35 + t10;
t5 = (-t109 + t35) * pkin(4) + t133;
t1 = -t108 * t10 + t106 * t9;
t2 = [0, 0, 0, 0, 0, qJDD(1), t140, t125, 0, 0, t96, 0.2e1 * t142, 0, t98, 0, 0, t105 * t72 + t173, t172 + (-t72 - t158) * t104, t139 + t174, pkin(1) * t72 + qJ(2) * t139, t96, 0, -0.2e1 * t142, 0, 0, t98, ((pkin(1) + 0.2e1 * t157 + 0.2e1 * t178) * qJDD(1) + t119) * t105 + t173, t105 * (pkin(2) * t87 + t175 + t71) + (qJ(3) * t87 + qJDD(3) + (qJD(1) * t85 + t136) * t104) * t104 + t174, -t172 + (t119 + 0.2e1 * t120) * t104, qJ(2) * (t104 * t33 + t105 * t129) + t126 * t40, t14, t104 * (-t106 * t52 - t162) + t105 * (-t108 * t52 + t168), t213, t117, -t216, t118, t104 * (qJ(3) * t50 + t171 - t209) + t105 * (t182 * t50 + t165 - t211) + t202, t104 * (t165 + t219) + t105 * (-t171 + t218) - t220 + t123 * t52, t104 * (-t3 - t208) + t105 * (-t4 - t210) + t201, t123 * t29 + t193 * (t104 * t3 + t105 * t4), t14, t213, t104 * (t106 * t32 + t162) + t105 * (t108 * t32 - t168), t118, t216, t117, t104 * (-t106 * t8 - t209) + t105 * (-t108 * t8 - t211) + (t104 * t137 + t105 * t135) * t50 + t202, t104 * (-t106 * t5 + t108 * t6 - t208) + t105 * (-t106 * t6 - t108 * t5 - t210) + t201, t104 * (t108 * t7 - t219) + t105 * (-t106 * t7 - t218) + t220 + (t104 * (-qJ(3) - t180) + t105 * (-t179 - t182) - pkin(1)) * t32, (t104 * (t137 + t180) + t105 * (t135 + t179) + pkin(1)) * t11 + t193 * (t104 * t1 + t105 * (t106 * t10 + t108 * t9)); 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t147, t97, -t87, -t72, 0, 0, 0, 0, 0, 0, -t147, -t87, -t97, -t40, 0, 0, 0, 0, 0, 0, -t50, -t52, t35, -t29, 0, 0, 0, 0, 0, 0, -t50, t35, t32, -t11; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t104 * t156, t97, -t153, t33, 0, 0, 0, 0, 0, 0, t196, -t16, t197, t3, 0, 0, 0, 0, 0, 0, t196, t197, t16, t1; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, t176, t190, t77, -t176, -t187, qJDD(4), -t12, -t13, 0, 0, t176, t121 + t159, -t190, qJDD(4), t187, -t176, pkin(4) * t204 + qJ(5) * t188 - t10, -pkin(4) * t77 - qJ(5) * t187, qJ(5) * t203 + (-t109 - t191) * pkin(4) + t133, -pkin(4) * t10 + qJ(5) * t9; 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, -t204, t77, t191, t10;];
tauJ_reg = t2;
