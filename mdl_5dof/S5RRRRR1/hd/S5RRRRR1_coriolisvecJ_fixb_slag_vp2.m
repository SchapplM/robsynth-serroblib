% Calculate vector of centrifugal and coriolis load on the joints for
% S5RRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 14:53
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S5RRRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 14:51:11
% EndTime: 2018-11-16 14:51:18
% DurationCPUTime: 3.40s
% Computational Cost: add. (6080->364), mult. (15877->539), div. (0->0), fcn. (12390->8), ass. (0->167)
t144 = qJD(2) + qJD(3);
t142 = qJD(4) + t144;
t146 = sin(qJ(5));
t150 = cos(qJ(5));
t148 = sin(qJ(3));
t149 = sin(qJ(2));
t152 = cos(qJ(3));
t153 = cos(qJ(2));
t131 = t148 * t149 - t152 * t153;
t129 = t131 * qJD(1);
t167 = t148 * t153 + t152 * t149;
t130 = t167 * qJD(1);
t147 = sin(qJ(4));
t151 = cos(qJ(4));
t171 = t129 * t147 - t151 * t130;
t74 = t142 * t150 - t146 * t171;
t241 = -t74 / 0.2e1;
t75 = t142 * t146 + t150 * t171;
t240 = -t75 / 0.2e1;
t183 = t151 * t129 + t130 * t147;
t87 = qJD(5) - t183;
t238 = -t87 / 0.2e1;
t105 = t144 * t131;
t98 = t105 * qJD(1);
t106 = t144 * t167;
t99 = t106 * qJD(1);
t42 = t183 * qJD(4) + t147 * t99 + t151 * t98;
t22 = qJD(5) * t74 + t150 * t42;
t23 = -qJD(5) * t75 - t146 * t42;
t8 = -mrSges(6,1) * t23 + mrSges(6,2) * t22;
t208 = pkin(2) * qJD(2);
t187 = t152 * t208;
t133 = pkin(3) * t144 + t187;
t196 = t148 * t151;
t169 = t147 * t152 + t196;
t193 = qJD(4) * t151;
t159 = (qJD(3) * t169 + t148 * t193) * pkin(2);
t194 = qJD(4) * t147;
t80 = qJD(2) * t159 + t133 * t194;
t256 = m(6) * t80 + t8;
t43 = qJD(4) * t171 + t147 * t98 - t151 * t99;
t195 = qJD(2) * t149;
t188 = pkin(2) * t195;
t81 = -pkin(3) * t99 - qJD(1) * t188;
t11 = pkin(4) * t43 - pkin(6) * t42 + t81;
t189 = t148 * t208;
t111 = t133 * t147 + t151 * t189;
t108 = pkin(6) * t142 + t111;
t209 = pkin(2) * qJD(1);
t134 = qJD(1) * pkin(1) + t153 * t209;
t109 = -pkin(3) * t129 + t134;
t46 = -pkin(4) * t183 - pkin(6) * t171 + t109;
t30 = -t108 * t146 + t150 * t46;
t197 = t147 * t148;
t168 = t151 * t152 - t197;
t160 = (qJD(3) * t168 - t148 * t194) * pkin(2);
t79 = qJD(2) * t160 + t133 * t193;
t2 = qJD(5) * t30 + t11 * t146 + t150 * t79;
t31 = t108 * t150 + t146 * t46;
t3 = -qJD(5) * t31 + t11 * t150 - t146 * t79;
t255 = -t3 * t146 + t150 * t2;
t110 = t133 * t151 - t147 * t189;
t107 = -pkin(4) * t142 - t110;
t179 = mrSges(6,1) * t146 + mrSges(6,2) * t150;
t165 = t107 * t179;
t172 = t146 * t31 + t150 * t30;
t174 = Ifges(6,5) * t150 - Ifges(6,6) * t146;
t212 = Ifges(6,4) * t150;
t176 = -Ifges(6,2) * t146 + t212;
t213 = Ifges(6,4) * t146;
t178 = Ifges(6,1) * t150 - t213;
t230 = t150 / 0.2e1;
t232 = -t146 / 0.2e1;
t239 = t75 / 0.2e1;
t226 = Ifges(6,4) * t75;
t28 = Ifges(6,2) * t74 + Ifges(6,6) * t87 + t226;
t73 = Ifges(6,4) * t74;
t29 = Ifges(6,1) * t75 + Ifges(6,5) * t87 + t73;
t254 = -t172 * mrSges(6,3) + t87 * t174 / 0.2e1 + t178 * t239 + t74 * t176 / 0.2e1 + t165 + t29 * t230 + t28 * t232;
t252 = (m(4) * t134 - mrSges(4,1) * t129 - mrSges(4,2) * t130) * pkin(2);
t227 = Ifges(5,4) * t171;
t86 = Ifges(5,4) * t183;
t251 = -t146 * t30 + t150 * t31;
t250 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t22 + Ifges(6,6) * t23;
t249 = Ifges(5,6) * t142 / 0.2e1 + t227 / 0.2e1;
t63 = pkin(4) * t171 - pkin(6) * t183;
t244 = Ifges(6,5) * t240 + Ifges(6,6) * t241 + Ifges(6,3) * t238;
t207 = t110 * mrSges(5,3);
t211 = Ifges(5,5) * t142;
t228 = Ifges(5,1) * t171;
t57 = t211 + t86 + t228;
t248 = t109 * mrSges(5,2) + t57 / 0.2e1 + t211 / 0.2e1 - t207 + t86 / 0.2e1 + t254;
t247 = -0.2e1 * pkin(1);
t246 = t22 / 0.2e1;
t245 = t23 / 0.2e1;
t243 = t43 / 0.2e1;
t224 = Ifges(5,2) * t183;
t242 = t224 / 0.2e1 + t249;
t235 = t129 / 0.2e1;
t234 = -t130 / 0.2e1;
t233 = t130 / 0.2e1;
t231 = t149 / 0.2e1;
t229 = t153 / 0.2e1;
t221 = pkin(3) * t130;
t140 = t153 * pkin(2) + pkin(1);
t218 = -mrSges(5,1) * t142 - mrSges(6,1) * t74 + mrSges(6,2) * t75 + mrSges(5,3) * t171;
t217 = mrSges(4,3) * t129;
t216 = mrSges(4,3) * t130;
t215 = Ifges(3,4) * t149;
t214 = Ifges(3,4) * t153;
t206 = t111 * t171;
t205 = t130 * Ifges(4,4);
t203 = t146 * t183;
t201 = t150 * t183;
t200 = Ifges(3,5) * qJD(2);
t199 = Ifges(3,6) * qJD(2);
t139 = pkin(2) * t152 + pkin(3);
t125 = pkin(2) * t196 + t147 * t139;
t192 = qJD(5) * t146;
t191 = qJD(5) * t150;
t190 = t149 * t209;
t116 = -pkin(3) * t131 + t140;
t181 = -t146 * t2 - t150 * t3;
t180 = mrSges(6,1) * t150 - mrSges(6,2) * t146;
t177 = Ifges(6,1) * t146 + t212;
t175 = Ifges(6,2) * t150 + t213;
t173 = Ifges(6,5) * t146 + Ifges(6,6) * t150;
t170 = t151 * t131 + t147 * t167;
t104 = t131 * t147 - t151 * t167;
t93 = -pkin(3) * t106 - t188;
t52 = -mrSges(6,2) * t87 + mrSges(6,3) * t74;
t53 = mrSges(6,1) * t87 - mrSges(6,3) * t75;
t77 = -mrSges(5,2) * t142 + mrSges(5,3) * t183;
t166 = -t146 * t53 + t150 * t52 + t77;
t124 = -pkin(2) * t197 + t139 * t151;
t51 = -t221 + t63;
t161 = m(6) * t172 + t146 * t52 + t150 * t53;
t10 = -mrSges(6,2) * t43 + mrSges(6,3) * t23;
t9 = mrSges(6,1) * t43 - mrSges(6,3) * t22;
t158 = -t53 * t191 - t52 * t192 - t146 * t9 + t150 * t10 + m(6) * (-t30 * t191 - t31 * t192 + t255);
t157 = -t109 * mrSges(5,1) - t30 * mrSges(6,1) + t31 * mrSges(6,2) + t111 * mrSges(5,3) + t242 + 0.2e1 * t244 + t249;
t6 = t22 * Ifges(6,4) + t23 * Ifges(6,2) + t43 * Ifges(6,6);
t7 = t22 * Ifges(6,1) + t23 * Ifges(6,4) + t43 * Ifges(6,5);
t156 = t177 * t246 + t175 * t245 + t173 * t243 + t146 * t7 / 0.2e1 + Ifges(5,5) * t42 + t6 * t230 - Ifges(5,6) * t43 - t79 * mrSges(5,2) + (-t180 - mrSges(5,1)) * t80 + t255 * mrSges(6,3) + t254 * qJD(5);
t123 = Ifges(4,4) * t129;
t84 = t129 * Ifges(4,2) + t144 * Ifges(4,6) - t205;
t85 = -t130 * Ifges(4,1) + t144 * Ifges(4,5) + t123;
t154 = -t134 * (-mrSges(4,1) * t130 + mrSges(4,2) * t129) - t31 * (-mrSges(6,2) * t171 - mrSges(6,3) * t203) - t30 * (mrSges(6,1) * t171 - mrSges(6,3) * t201) + t171 * t242 + t171 * t244 + (Ifges(6,3) * t171 + t174 * t183) * t238 + (Ifges(6,5) * t171 + t178 * t183) * t240 + (Ifges(6,6) * t171 + t176 * t183) * t241 - t142 * (Ifges(5,5) * t183 - Ifges(5,6) * t171) / 0.2e1 + t183 * t207 - t183 * t165 + (Ifges(4,1) * t129 + t205) * t233 + t84 * t234 + t187 * t217 - t144 * (Ifges(4,5) * t129 + Ifges(4,6) * t130) / 0.2e1 + Ifges(4,5) * t98 + Ifges(4,6) * t99 - t171 * (Ifges(5,1) * t183 - t227) / 0.2e1 + t156 - t29 * t201 / 0.2e1 + t28 * t203 / 0.2e1 - (-Ifges(5,2) * t171 + t57 + t86) * t183 / 0.2e1 - (Ifges(4,2) * t130 + t123 + t85) * t129 / 0.2e1;
t128 = t200 + (-Ifges(3,1) * t149 - t214) * qJD(1);
t127 = t199 + (-Ifges(3,2) * t153 - t215) * qJD(1);
t122 = t168 * t208;
t121 = t169 * t208;
t119 = -pkin(4) - t124;
t114 = mrSges(4,1) * t144 + t216;
t113 = -mrSges(4,2) * t144 + t217;
t112 = -t190 - t221;
t97 = t139 * t194 + t159;
t62 = -mrSges(5,1) * t183 + mrSges(5,2) * t171;
t61 = mrSges(5,1) * t171 + mrSges(5,2) * t183;
t50 = qJD(4) * t104 + t105 * t147 - t151 * t106;
t49 = qJD(4) * t170 + t105 * t151 + t106 * t147;
t45 = t110 * t150 + t146 * t63;
t44 = -t110 * t146 + t150 * t63;
t39 = Ifges(6,3) * t43;
t36 = t122 * t150 + t146 * t51;
t35 = -t122 * t146 + t150 * t51;
t1 = [t140 * (-mrSges(4,1) * t99 + mrSges(4,2) * t98) + t144 * (Ifges(4,5) * t105 + Ifges(4,6) * t106) / 0.2e1 + t134 * (-mrSges(4,1) * t106 + mrSges(4,2) * t105) + t116 * (mrSges(5,1) * t43 + mrSges(5,2) * t42) + t105 * t85 / 0.2e1 + t106 * t84 / 0.2e1 + t93 * t62 + (t106 * t235 + t131 * t99) * Ifges(4,2) + (t105 * t234 - t167 * t98) * Ifges(4,1) + (t109 * t93 + t116 * t81) * m(5) + (-t224 / 0.2e1 - t157) * t50 + (t228 / 0.2e1 + t248) * t49 + t161 * (pkin(4) * t50 - pkin(6) * t49 + t93) + (t52 * t191 - t53 * t192 + t146 * t10 + t150 * t9 + m(6) * (t31 * t191 - t30 * t192 - t181)) * (-pkin(4) * t170 - pkin(6) * t104 + t116) + ((-t128 / 0.2e1 - t200 / 0.2e1 + (mrSges(3,2) * t247 + 0.3e1 / 0.2e1 * t214) * qJD(1)) * t153 + (-t105 * t152 + t106 * t148 + (t131 * t152 - t148 * t167) * qJD(3)) * pkin(2) * mrSges(4,3)) * qJD(2) + (t105 * t235 + t106 * t234 + t131 * t98 - t167 * t99) * Ifges(4,4) - (t81 * mrSges(5,1) - Ifges(5,4) * t42 + t39 / 0.2e1 - t79 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2)) * t43 + t250) * t170 + (t174 * t243 + t176 * t245 + t178 * t246 + t81 * mrSges(5,2) + Ifges(5,1) * t42 - Ifges(5,4) * t43 + t6 * t232 + t7 * t230 + (mrSges(5,3) + t179) * t80 + t181 * mrSges(6,3) + (-t150 * t28 / 0.2e1 + t29 * t232 + t173 * t238 + t175 * t241 + t177 * t240 + t107 * t180 - t251 * mrSges(6,3)) * qJD(5)) * t104 + (t127 / 0.2e1 + t199 / 0.2e1 - t252 + (mrSges(3,1) * t247 - 0.3e1 / 0.2e1 * t215 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t153 + (-m(4) * t140 + mrSges(4,1) * t131 + mrSges(4,2) * t167) * pkin(2)) * qJD(1)) * t195; (-t149 * t127 / 0.2e1 + t128 * t229 + ((-Ifges(3,1) * t153 + t215) * t231 + (Ifges(3,2) * t149 - t214) * t229 - pkin(1) * (-mrSges(3,1) * t149 - mrSges(3,2) * t153)) * qJD(1) + t149 * t252 + (-Ifges(3,5) * t153 / 0.2e1 + Ifges(3,6) * t231) * qJD(2)) * qJD(1) + t218 * t97 + t154 + ((-t152 * t98 + (-qJD(2) * t130 + t99) * t148) * mrSges(4,3) + (t113 * t152 - t114 * t148 + (-mrSges(4,1) * t148 - mrSges(4,2) * t152) * qJD(2)) * qJD(3)) * pkin(2) - t112 * t62 + t119 * t8 + m(6) * (t107 * t97 + t119 * t80) + m(5) * (-t110 * t97 - t124 * t80 + t125 * t79) + (-m(5) * t112 - t61) * t109 + (-t124 * t42 - t125 * t43 + t206) * mrSges(5,3) - t161 * (t51 - t190) + t158 * (pkin(6) + t125) + (m(5) * t111 + m(6) * t251 + t166) * (t139 * t193 + t160); -t36 * t52 - t35 * t53 + t154 - m(6) * (t107 * t121 + t30 * t35 + t31 * t36) - m(5) * (-t110 * t121 + t111 * t122) - t122 * t77 - t109 * t61 - t218 * t121 + ((-mrSges(4,2) * qJD(3) - t113) * t152 + (-mrSges(4,1) * qJD(3) + t114 - t216) * t148) * t208 + mrSges(5,3) * t206 + (t130 * t62 + (-t147 * t43 - t151 * t42) * mrSges(5,3) + (t218 * t147 + t166 * t151 + m(6) * (t107 * t147 + t251 * t151)) * qJD(4) + (0.2e1 * t109 * t233 + t147 * t79 - t151 * t80 + (-t110 * t147 + t111 * t151) * qJD(4)) * m(5)) * pkin(3) + t158 * (pkin(3) * t147 + pkin(6)) + t256 * (-pkin(3) * t151 - pkin(4)); -t218 * t111 - t45 * t52 - t44 * t53 + t157 * t171 - t110 * t77 + ((-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t171 - t248) * t183 - m(6) * (t107 * t111 + t30 * t44 + t31 * t45) + t156 + t158 * pkin(6) - t256 * pkin(4); t39 - t107 * (mrSges(6,1) * t75 + mrSges(6,2) * t74) + (Ifges(6,1) * t74 - t226) * t240 + t28 * t239 + (Ifges(6,5) * t74 - Ifges(6,6) * t75) * t238 - t30 * t52 + t31 * t53 + (t30 * t74 + t31 * t75) * mrSges(6,3) + (-Ifges(6,2) * t75 + t29 + t73) * t241 + t250;];
tauc  = t1(:);
