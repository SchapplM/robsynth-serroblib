% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 19:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR10_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR10_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:09:43
% EndTime: 2019-12-31 19:09:57
% DurationCPUTime: 5.42s
% Computational Cost: add. (6245->415), mult. (16359->587), div. (0->0), fcn. (12103->8), ass. (0->192)
t160 = sin(qJ(4));
t239 = -pkin(8) - pkin(7);
t191 = qJD(4) * t239;
t157 = sin(pkin(9));
t158 = cos(pkin(9));
t161 = sin(qJ(3));
t164 = cos(qJ(3));
t267 = -t157 * t161 + t164 * t158;
t134 = t267 * qJD(1);
t197 = t160 * t134;
t140 = t157 * t164 + t158 * t161;
t135 = t140 * qJD(1);
t107 = pkin(3) * t135 - pkin(7) * t134;
t223 = pkin(6) + qJ(2);
t148 = t223 * t157;
t141 = qJD(1) * t148;
t149 = t223 * t158;
t142 = qJD(1) * t149;
t110 = -t141 * t164 - t161 * t142;
t163 = cos(qJ(4));
t54 = t160 * t107 + t163 * t110;
t272 = pkin(8) * t197 + t160 * t191 - t54;
t226 = pkin(8) * t163;
t53 = t163 * t107 - t110 * t160;
t271 = -pkin(4) * t135 + t134 * t226 + t163 * t191 - t53;
t105 = -qJD(3) * pkin(3) - t110;
t121 = qJD(3) * t163 - t135 * t160;
t132 = qJD(4) - t134;
t111 = -t141 * t161 + t142 * t164;
t106 = qJD(3) * pkin(7) + t111;
t192 = -pkin(2) * t158 - pkin(1);
t147 = qJD(1) * t192 + qJD(2);
t84 = -pkin(3) * t134 - pkin(7) * t135 + t147;
t47 = -t106 * t160 + t163 * t84;
t48 = t106 * t163 + t160 * t84;
t177 = t160 * t48 + t163 * t47;
t180 = Ifges(5,5) * t163 - Ifges(5,6) * t160;
t220 = Ifges(5,4) * t163;
t182 = -Ifges(5,2) * t160 + t220;
t221 = Ifges(5,4) * t160;
t184 = Ifges(5,1) * t163 - t221;
t185 = mrSges(5,1) * t160 + mrSges(5,2) * t163;
t228 = t163 / 0.2e1;
t229 = -t160 / 0.2e1;
t122 = qJD(3) * t160 + t135 * t163;
t234 = t122 / 0.2e1;
t222 = Ifges(5,4) * t122;
t56 = t121 * Ifges(5,2) + t132 * Ifges(5,6) + t222;
t118 = Ifges(5,4) * t121;
t57 = t122 * Ifges(5,1) + t132 * Ifges(5,5) + t118;
t270 = t132 * t180 / 0.2e1 + t184 * t234 + t121 * t182 / 0.2e1 + t105 * t185 + t56 * t229 + t57 * t228 - t177 * mrSges(5,3);
t162 = cos(qJ(5));
t159 = sin(qJ(5));
t41 = pkin(8) * t121 + t48;
t211 = t159 * t41;
t40 = -pkin(8) * t122 + t47;
t34 = pkin(4) * t132 + t40;
t10 = t162 * t34 - t211;
t209 = t162 * t41;
t11 = t159 * t34 + t209;
t137 = t140 * qJD(3);
t130 = qJD(1) * t137;
t126 = Ifges(6,3) * t130;
t188 = t162 * t121 - t122 * t159;
t70 = t121 * t159 + t122 * t162;
t227 = Ifges(6,4) * t70;
t128 = qJD(5) + t132;
t233 = -t128 / 0.2e1;
t243 = -t70 / 0.2e1;
t245 = -t188 / 0.2e1;
t65 = Ifges(6,4) * t188;
t33 = Ifges(6,1) * t70 + Ifges(6,5) * t128 + t65;
t62 = -pkin(4) * t121 + t105;
t269 = t126 + (Ifges(6,5) * t188 - Ifges(6,6) * t70) * t233 + (t10 * t188 + t11 * t70) * mrSges(6,3) + (-Ifges(6,2) * t70 + t33 + t65) * t245 - t62 * (mrSges(6,1) * t70 + mrSges(6,2) * t188) + (Ifges(6,1) * t188 - t227) * t243;
t131 = Ifges(4,4) * t134;
t252 = t131 / 0.2e1 + t135 * Ifges(4,1) / 0.2e1;
t268 = t147 * mrSges(4,2) + Ifges(4,5) * qJD(3) + t252 + t270;
t193 = qJD(4) * t163;
t194 = qJD(4) * t160;
t169 = t267 * qJD(2);
t71 = qJD(1) * t169 + qJD(3) * t110;
t136 = t267 * qJD(3);
t129 = qJD(1) * t136;
t92 = pkin(3) * t130 - pkin(7) * t129;
t21 = -t106 * t194 + t160 * t92 + t163 * t71 + t84 * t193;
t79 = -qJD(4) * t122 - t129 * t160;
t12 = pkin(8) * t79 + t21;
t22 = -qJD(4) * t48 - t160 * t71 + t163 * t92;
t78 = qJD(4) * t121 + t129 * t163;
t9 = pkin(4) * t130 - pkin(8) * t78 + t22;
t2 = qJD(5) * t10 + t12 * t162 + t159 * t9;
t27 = qJD(5) * t188 + t159 * t79 + t162 * t78;
t28 = -qJD(5) * t70 - t159 * t78 + t162 * t79;
t3 = -qJD(5) * t11 - t12 * t159 + t162 * t9;
t266 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t27 + Ifges(6,6) * t28;
t247 = t27 / 0.2e1;
t246 = t28 / 0.2e1;
t32 = Ifges(6,2) * t188 + Ifges(6,6) * t128 + t227;
t264 = t32 / 0.2e1;
t231 = t130 / 0.2e1;
t150 = t239 * t160;
t151 = t239 * t163;
t119 = t150 * t162 + t151 * t159;
t263 = qJD(5) * t119 + t271 * t159 + t272 * t162;
t120 = t150 * t159 - t151 * t162;
t262 = -qJD(5) * t120 - t272 * t159 + t271 * t162;
t174 = t159 * t160 - t162 * t163;
t254 = qJD(4) + qJD(5);
t113 = t254 * t174;
t91 = t174 * t134;
t207 = -t113 + t91;
t144 = t159 * t163 + t160 * t162;
t114 = t254 * t144;
t90 = t144 * t134;
t206 = -t114 + t90;
t39 = -mrSges(6,1) * t188 + mrSges(6,2) * t70;
t257 = m(6) * t62 + t39;
t101 = t174 * t140;
t109 = -pkin(3) * t267 - pkin(7) * t140 + t192;
t117 = -t148 * t161 + t149 * t164;
t112 = t163 * t117;
t61 = t160 * t109 + t112;
t256 = -t164 * t148 - t149 * t161;
t255 = -t160 * t22 + t163 * t21;
t251 = (m(3) * qJ(2) + mrSges(3,3)) * (t157 ^ 2 + t158 ^ 2);
t250 = t22 * mrSges(5,1) - t21 * mrSges(5,2) + Ifges(5,5) * t78 + Ifges(5,6) * t79 + t266;
t249 = Ifges(6,4) * t247 + Ifges(6,2) * t246 + Ifges(6,6) * t231;
t248 = Ifges(6,1) * t247 + Ifges(6,4) * t246 + Ifges(6,5) * t231;
t244 = t188 / 0.2e1;
t242 = t70 / 0.2e1;
t241 = t78 / 0.2e1;
t240 = t79 / 0.2e1;
t236 = -t121 / 0.2e1;
t235 = -t122 / 0.2e1;
t232 = t128 / 0.2e1;
t230 = -t132 / 0.2e1;
t170 = t140 * qJD(2);
t72 = qJD(1) * t170 + qJD(3) * t111;
t219 = t256 * t72;
t214 = t134 * Ifges(4,2);
t205 = qJD(3) * mrSges(4,1) + mrSges(5,1) * t121 - mrSges(5,2) * t122 - t135 * mrSges(4,3);
t201 = t140 * t160;
t108 = pkin(3) * t137 - pkin(7) * t136;
t85 = qJD(3) * t256 + t169;
t190 = t163 * t108 - t160 * t85;
t189 = t130 * mrSges(4,1) + t129 * mrSges(4,2);
t60 = t163 * t109 - t117 * t160;
t186 = mrSges(5,1) * t163 - mrSges(5,2) * t160;
t183 = Ifges(5,1) * t160 + t220;
t181 = Ifges(5,2) * t163 + t221;
t179 = Ifges(5,5) * t160 + Ifges(5,6) * t163;
t45 = -pkin(4) * t267 - t140 * t226 + t60;
t49 = -pkin(8) * t201 + t61;
t19 = -t159 * t49 + t162 * t45;
t20 = t159 * t45 + t162 * t49;
t178 = -t160 * t21 - t163 * t22;
t176 = t160 * t47 - t163 * t48;
t82 = -mrSges(5,2) * t132 + mrSges(5,3) * t121;
t83 = mrSges(5,1) * t132 - mrSges(5,3) * t122;
t175 = -t160 * t83 + t163 * t82;
t172 = t136 * t160 + t140 * t193;
t29 = t160 * t108 + t109 * t193 - t117 * t194 + t163 * t85;
t86 = qJD(3) * t117 + t170;
t167 = t11 * mrSges(6,2) + t48 * mrSges(5,2) - t128 * Ifges(6,3) - t70 * Ifges(6,5) - t188 * Ifges(6,6) - t132 * Ifges(5,3) - t122 * Ifges(5,5) - t121 * Ifges(5,6) + Ifges(4,6) * qJD(3) + t135 * Ifges(4,4) + t214 / 0.2e1 - t10 * mrSges(6,1) - t147 * mrSges(4,1) - t47 * mrSges(5,1);
t154 = -pkin(4) * t163 - pkin(3);
t127 = Ifges(5,3) * t130;
t123 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t134;
t100 = t144 * t140;
t87 = pkin(4) * t201 - t256;
t74 = pkin(4) * t197 + t111;
t59 = -mrSges(5,2) * t130 + mrSges(5,3) * t79;
t58 = mrSges(5,1) * t130 - mrSges(5,3) * t78;
t52 = mrSges(6,1) * t128 - mrSges(6,3) * t70;
t51 = -mrSges(6,2) * t128 + mrSges(6,3) * t188;
t50 = pkin(4) * t172 + t86;
t44 = -mrSges(5,1) * t79 + mrSges(5,2) * t78;
t43 = -pkin(4) * t79 + t72;
t38 = t78 * Ifges(5,1) + t79 * Ifges(5,4) + t130 * Ifges(5,5);
t37 = t78 * Ifges(5,4) + t79 * Ifges(5,2) + t130 * Ifges(5,6);
t36 = t101 * t254 - t144 * t136;
t35 = -t114 * t140 - t136 * t174;
t30 = -qJD(4) * t61 + t190;
t24 = -mrSges(6,2) * t130 + mrSges(6,3) * t28;
t23 = mrSges(6,1) * t130 - mrSges(6,3) * t27;
t18 = -pkin(8) * t172 + t29;
t15 = -t136 * t226 + pkin(4) * t137 + (-t112 + (pkin(8) * t140 - t109) * t160) * qJD(4) + t190;
t14 = t162 * t40 - t211;
t13 = -t159 * t40 - t209;
t8 = -mrSges(6,1) * t28 + mrSges(6,2) * t27;
t5 = -qJD(5) * t20 + t15 * t162 - t159 * t18;
t4 = qJD(5) * t19 + t15 * t159 + t162 * t18;
t1 = [(-Ifges(6,4) * t101 - Ifges(6,2) * t100) * t246 + (-Ifges(6,1) * t101 - Ifges(6,4) * t100) * t247 + (-Ifges(6,5) * t101 - Ifges(6,6) * t100) * t231 + (-t10 * t35 - t100 * t2 + t101 * t3 + t11 * t36) * mrSges(6,3) + t43 * (mrSges(6,1) * t100 - mrSges(6,2) * t101) + t36 * t264 - t256 * t44 - t205 * t86 + (-t214 / 0.2e1 - t167) * t137 + m(5) * (t105 * t86 + t21 * t61 + t22 * t60 + t29 * t48 + t30 * t47 - t219) + m(4) * (-t110 * t86 + t111 * t85 + t117 * t71 - t219) + 0.2e1 * t251 * qJD(2) * qJD(1) + (Ifges(6,1) * t35 + Ifges(6,4) * t36) * t242 + (Ifges(6,4) * t35 + Ifges(6,2) * t36) * t244 - t101 * t248 - t100 * t249 + (Ifges(6,5) * t35 + Ifges(6,6) * t36) * t232 + (t252 + t268) * t136 + t85 * t123 + t29 * t82 + t30 * t83 + t87 * t8 + t62 * (-mrSges(6,1) * t36 + mrSges(6,2) * t35) + t60 * t58 + t61 * t59 + t50 * t39 + t4 * t51 + t5 * t52 + m(6) * (t10 * t5 + t11 * t4 + t19 * t3 + t2 * t20 + t43 * t87 + t50 * t62) + t35 * t33 / 0.2e1 + t19 * t23 + t20 * t24 + (t180 * t231 + t184 * t241 + t182 * t240 - Ifges(4,4) * t130 + Ifges(4,1) * t129 + t37 * t229 + t38 * t228 + (mrSges(4,3) + t185) * t72 + t178 * mrSges(5,3) + (t105 * t186 + t179 * t230 + t183 * t235 + t181 * t236 + t57 * t229 - t163 * t56 / 0.2e1 + t176 * mrSges(5,3)) * qJD(4)) * t140 - (-Ifges(4,4) * t129 + t126 / 0.2e1 + t127 / 0.2e1 - t71 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) + Ifges(6,3) / 0.2e1) * t130 + t250) * t267 + (-t110 * t136 - t111 * t137 - t117 * t130 - t129 * t256) * mrSges(4,3) + t192 * t189; -t174 * t23 + t144 * t24 + t160 * t59 + t163 * t58 + t206 * t52 + t207 * t51 + t175 * qJD(4) + (-t39 + t205) * t135 + (-t123 - t175) * t134 - m(4) * (-t110 * t135 + t111 * t134) + t189 - t251 * qJD(1) ^ 2 + (t10 * t206 + t11 * t207 - t135 * t62 + t144 * t2 - t174 * t3) * m(6) + (-t105 * t135 - t132 * t176 - t178) * m(5); (-mrSges(4,1) - t186) * t72 + (-Ifges(6,1) * t91 - Ifges(6,4) * t90) * t243 + (-Ifges(6,4) * t91 - Ifges(6,2) * t90) * t245 + (-Ifges(6,5) * t91 - Ifges(6,6) * t90) * t233 + (t111 * mrSges(4,3) + t167) * t135 + (m(5) * t255 - t160 * t58 + t163 * t59 + (-m(5) * t177 - t160 * t82 - t163 * t83) * qJD(4)) * pkin(7) + t255 * mrSges(5,3) + t205 * t111 + (-pkin(3) * t72 - t105 * t111 - t47 * t53 - t48 * t54) * m(5) + (pkin(4) * t160 * t257 + t270) * qJD(4) + (Ifges(6,5) * t144 - Ifges(6,6) * t174 + t179) * t231 + (Ifges(6,4) * t144 - Ifges(6,2) * t174) * t246 + (Ifges(6,1) * t144 - Ifges(6,4) * t174) * t247 + t43 * (mrSges(6,1) * t174 + mrSges(6,2) * t144) + (-t10 * t207 + t11 * t206 - t144 * t3 - t174 * t2) * mrSges(6,3) - t174 * t249 + (-t113 / 0.2e1 + t91 / 0.2e1) * t33 + (-t114 / 0.2e1 + t90 / 0.2e1) * t32 + (-Ifges(6,1) * t113 - Ifges(6,4) * t114) * t242 + (-Ifges(6,4) * t113 - Ifges(6,2) * t114) * t244 + (-Ifges(6,5) * t113 - Ifges(6,6) * t114) * t232 + t181 * t240 + t183 * t241 + t144 * t248 + t37 * t228 + t262 * t52 + (t10 * t262 + t11 * t263 + t119 * t3 + t120 * t2 + t154 * t43 - t62 * t74) * m(6) + t263 * t51 + t160 * t38 / 0.2e1 + t154 * t8 + (-t131 / 0.2e1 + t110 * mrSges(4,3) + (-Ifges(4,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t135 - t268) * t134 - Ifges(4,6) * t130 + Ifges(4,5) * t129 - t110 * t123 + t119 * t23 + t120 * t24 - t54 * t82 - t53 * t83 - t71 * mrSges(4,2) - t74 * t39 - pkin(3) * t44 + (-mrSges(6,1) * t206 + mrSges(6,2) * t207) * t62; t127 - m(6) * (t10 * t13 + t11 * t14) + (t159 * t24 + t162 * t23 + m(6) * (t159 * t2 + t162 * t3) - t257 * t122 + (-t159 * t52 + t162 * t51 + m(6) * (-t10 * t159 + t11 * t162)) * qJD(5)) * pkin(4) + (-Ifges(5,2) * t122 + t118 + t57) * t236 + t250 + t70 * t264 + (Ifges(5,5) * t121 - Ifges(5,6) * t122) * t230 + t56 * t234 + (Ifges(5,1) * t121 - t222) * t235 - t105 * (mrSges(5,1) * t122 + mrSges(5,2) * t121) - t47 * t82 + t48 * t83 - t14 * t51 - t13 * t52 + (t121 * t47 + t122 * t48) * mrSges(5,3) + t269; -t10 * t51 + t11 * t52 + t32 * t242 + t266 + t269;];
tauc = t1(:);
