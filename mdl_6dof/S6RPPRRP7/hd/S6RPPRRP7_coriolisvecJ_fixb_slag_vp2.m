% Calculate vector of centrifugal and coriolis load on the joints for
% S6RPPRRP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [6x1]
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-23 15:47
% Revision: 76f9d5e39f14dc242b53c0d9d3d9db48bd8f37c0
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S6RPPRRP7_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert( isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP7_coriolisvecJ_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-23 15:46:53
% EndTime: 2018-11-23 15:46:58
% DurationCPUTime: 4.71s
% Computational Cost: add. (4768->411), mult. (10964->523), div. (0->0), fcn. (7379->6), ass. (0->190)
t269 = Ifges(6,4) + Ifges(7,4);
t272 = -Ifges(5,1) / 0.2e1;
t270 = Ifges(6,1) + Ifges(7,1);
t262 = Ifges(6,5) + Ifges(7,5);
t268 = Ifges(6,2) + Ifges(7,2);
t261 = Ifges(7,6) + Ifges(6,6);
t137 = sin(pkin(9));
t138 = cos(pkin(9));
t231 = sin(qJ(4));
t232 = cos(qJ(4));
t112 = t137 * t231 - t138 * t232;
t271 = Ifges(5,2) / 0.2e1;
t108 = t112 * qJD(1);
t140 = sin(qJ(5));
t141 = cos(qJ(5));
t91 = qJD(4) * t141 + t108 * t140;
t267 = t269 * t91;
t92 = qJD(4) * t140 - t108 * t141;
t266 = t269 * t92;
t113 = t137 * t232 + t138 * t231;
t107 = t113 * qJD(1);
t136 = qJD(1) * qJ(2);
t132 = qJD(3) + t136;
t133 = t137 * pkin(3);
t117 = qJD(1) * t133 + t132;
t139 = -pkin(1) - qJ(3);
t122 = qJD(1) * t139 + qJD(2);
t179 = -pkin(7) * qJD(1) + t122;
t101 = t179 * t137;
t102 = t179 * t138;
t72 = t101 * t232 + t102 * t231;
t64 = qJD(4) * pkin(8) + t72;
t65 = pkin(4) * t107 + pkin(8) * t108 + t117;
t22 = -t140 * t64 + t141 * t65;
t23 = t140 * t65 + t141 * t64;
t151 = t140 * t23 + t22 * t141;
t153 = Ifges(7,5) * t141 - Ifges(7,6) * t140;
t155 = Ifges(6,5) * t141 - Ifges(6,6) * t140;
t212 = Ifges(7,4) * t141;
t157 = -Ifges(7,2) * t140 + t212;
t214 = Ifges(6,4) * t141;
t159 = -Ifges(6,2) * t140 + t214;
t213 = Ifges(7,4) * t140;
t161 = Ifges(7,1) * t141 - t213;
t215 = Ifges(6,4) * t140;
t163 = Ifges(6,1) * t141 - t215;
t164 = mrSges(7,1) * t140 + mrSges(7,2) * t141;
t166 = mrSges(6,1) * t140 + mrSges(6,2) * t141;
t11 = qJ(6) * t91 + t23;
t10 = -qJ(6) * t92 + t22;
t255 = qJD(5) + t107;
t9 = pkin(5) * t255 + t10;
t169 = t11 * t140 + t9 * t141;
t233 = t141 / 0.2e1;
t235 = t140 / 0.2e1;
t237 = t255 / 0.2e1;
t240 = t92 / 0.2e1;
t242 = t91 / 0.2e1;
t251 = t255 * t262 + t270 * t92 + t267;
t260 = t255 * t261 + t268 * t91 + t266;
t71 = -t101 * t231 + t102 * t232;
t63 = -qJD(4) * pkin(4) - t71;
t28 = -pkin(5) * t91 + qJD(6) + t63;
t248 = t151 * mrSges(6,3) + t169 * mrSges(7,3) - t164 * t28 - t166 * t63 - (t157 + t159) * t242 - (t161 + t163) * t240 - (t153 + t155) * t237 + t260 * t235 - t251 * t233;
t264 = t108 * t272;
t265 = t117 * mrSges(5,2) - t71 * mrSges(5,3) - Ifges(5,4) * t107 + Ifges(5,5) * qJD(4) - t248 + t264;
t263 = t107 * t271;
t259 = t112 * qJD(3);
t258 = t140 * t262 + t261 * t141;
t257 = t268 * t141 + t213 + t215;
t256 = t140 * t270 + t212 + t214;
t243 = -t91 / 0.2e1;
t241 = -t92 / 0.2e1;
t238 = -t255 / 0.2e1;
t254 = m(6) * t151;
t100 = t108 * qJD(4);
t99 = qJD(4) * t107;
t55 = qJD(5) * t91 - t141 * t99;
t56 = -qJD(5) * t92 + t140 * t99;
t253 = -t100 * t261 + t268 * t56 + t269 * t55;
t252 = -t100 * t262 + t269 * t56 + t270 * t55;
t198 = t137 ^ 2 + t138 ^ 2;
t250 = mrSges(4,3) * t198;
t223 = -pkin(7) + t139;
t115 = t223 * t137;
t116 = t223 * t138;
t87 = t115 * t231 - t116 * t232;
t196 = qJD(5) * t141;
t197 = qJD(5) * t140;
t146 = t113 * qJD(3);
t42 = -qJD(1) * t146 + qJD(4) * t71;
t195 = qJD(1) * qJD(2);
t70 = -pkin(4) * t100 + pkin(8) * t99 + t195;
t5 = t140 * t70 + t141 * t42 + t196 * t65 - t197 * t64;
t6 = -qJD(5) * t23 - t140 * t42 + t141 * t70;
t172 = -t140 * t6 + t141 * t5;
t190 = Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1;
t191 = -Ifges(7,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t192 = Ifges(7,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t246 = t190 * t255 - t191 * t91 + t192 * t92 - Ifges(5,6) * qJD(4) - t11 * mrSges(7,2) + t9 * mrSges(7,1) - t72 * mrSges(5,3) + t117 * mrSges(5,1) + t108 * Ifges(5,4) + t263 - t23 * mrSges(6,2) + t22 * mrSges(6,1) - t261 * t243 - t262 * t241 - (Ifges(7,3) + Ifges(6,3)) * t238;
t230 = m(3) * qJ(2);
t43 = -qJD(1) * t259 + qJD(4) * t72;
t227 = t43 * t87;
t222 = -qJ(6) - pkin(8);
t37 = -mrSges(7,1) * t100 - mrSges(7,3) * t55;
t38 = -mrSges(6,1) * t100 - mrSges(6,3) * t55;
t221 = -t37 - t38;
t39 = mrSges(7,2) * t100 + mrSges(7,3) * t56;
t40 = mrSges(6,2) * t100 + mrSges(6,3) * t56;
t220 = t39 + t40;
t219 = -qJD(4) * mrSges(5,1) - mrSges(6,1) * t91 + mrSges(6,2) * t92 - mrSges(5,3) * t108;
t85 = -pkin(4) * t108 + pkin(8) * t107;
t27 = t140 * t85 + t141 * t71;
t66 = -mrSges(7,2) * t255 + mrSges(7,3) * t91;
t67 = -mrSges(6,2) * t255 + mrSges(6,3) * t91;
t218 = t66 + t67;
t68 = mrSges(7,1) * t255 - mrSges(7,3) * t92;
t69 = mrSges(6,1) * t255 - mrSges(6,3) * t92;
t217 = t68 + t69;
t128 = qJ(2) + t133;
t84 = pkin(4) * t113 + pkin(8) * t112 + t128;
t88 = t115 * t232 + t116 * t231;
t86 = t141 * t88;
t30 = t140 * t84 + t86;
t216 = m(4) * qJD(3);
t208 = t112 * t43;
t168 = mrSges(4,1) * t137 + mrSges(4,2) * t138;
t206 = mrSges(5,1) * t107 - mrSges(5,2) * t108 + qJD(1) * t168;
t203 = qJ(6) * t141;
t202 = t107 * t140;
t183 = qJD(4) * t231;
t184 = qJD(4) * t232;
t110 = -t137 * t183 + t138 * t184;
t201 = t110 * t140;
t200 = t110 * t141;
t199 = t112 * t140;
t45 = -mrSges(7,1) * t91 + mrSges(7,2) * t92;
t194 = -t45 - t219;
t57 = -qJD(4) * t87 - t146;
t109 = -t137 * t184 - t138 * t183;
t82 = pkin(4) * t110 - pkin(8) * t109 + qJD(2);
t193 = t140 * t82 + t141 * t57 + t196 * t84;
t185 = t112 * t196;
t20 = -mrSges(7,1) * t56 + mrSges(7,2) * t55;
t182 = -t100 * mrSges(5,1) - mrSges(5,2) * t99;
t181 = -t140 * t57 + t141 * t82;
t26 = -t140 * t71 + t141 * t85;
t29 = -t140 * t88 + t141 * t84;
t180 = qJD(5) * t222;
t178 = qJD(1) * t198;
t1 = -pkin(5) * t100 - qJ(6) * t55 - qJD(6) * t92 + t6;
t2 = qJ(6) * t56 + qJD(6) * t91 + t5;
t174 = t1 * t141 + t140 * t2;
t173 = -t1 * t140 + t141 * t2;
t171 = t140 * t5 + t141 * t6;
t170 = t11 * t141 - t140 * t9;
t167 = -mrSges(6,1) * t141 + mrSges(6,2) * t140;
t165 = -mrSges(7,1) * t141 + mrSges(7,2) * t140;
t150 = -t140 * t22 + t141 * t23;
t149 = -qJ(6) * t109 + qJD(6) * t112;
t148 = -t140 * t218 - t141 * t217;
t147 = t6 * mrSges(6,1) + t1 * mrSges(7,1) - t5 * mrSges(6,2) - t2 * mrSges(7,2);
t58 = qJD(4) * t88 - t259;
t142 = qJD(1) ^ 2;
t131 = -pkin(5) * t141 - pkin(4);
t119 = t222 * t141;
t118 = t222 * t140;
t106 = -qJD(6) * t140 + t141 * t180;
t105 = qJD(6) * t141 + t140 * t180;
t98 = Ifges(6,3) * t100;
t97 = Ifges(7,3) * t100;
t94 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t107;
t61 = -pkin(5) * t199 + t87;
t52 = Ifges(6,5) * t55;
t51 = Ifges(7,5) * t55;
t50 = Ifges(6,6) * t56;
t49 = Ifges(7,6) * t56;
t44 = -pkin(5) * t202 + t72;
t25 = (t109 * t140 - t185) * pkin(5) + t58;
t24 = qJ(6) * t199 + t30;
t21 = -mrSges(6,1) * t56 + mrSges(6,2) * t55;
t19 = pkin(5) * t113 + t112 * t203 + t29;
t18 = qJ(6) * t202 + t27;
t17 = -pkin(5) * t56 + t43;
t12 = -pkin(5) * t108 + t107 * t203 + t26;
t8 = -qJD(5) * t30 + t181;
t7 = -t197 * t88 + t193;
t4 = qJ(6) * t185 + (-qJD(5) * t88 + t149) * t140 + t193;
t3 = pkin(5) * t110 + t149 * t141 + (-t86 + (-qJ(6) * t112 - t84) * t140) * qJD(5) + t181;
t13 = [(t263 + t246) * t110 + m(6) * (t22 * t8 + t23 * t7 + t29 * t6 + t30 * t5 + t58 * t63 + t227) + m(5) * (t42 * t88 + t57 * t72 - t58 * t71 + t227) + t219 * t58 + (t264 + t265) * t109 + (((2 * mrSges(3,3)) + t168 + 0.2e1 * t230) * qJD(1) + t206 + m(5) * (qJD(1) * t128 + t117) + m(4) * (t132 + t136)) * qJD(2) + m(7) * (t1 * t19 + t11 * t4 + t17 * t61 + t2 * t24 + t25 * t28 + t3 * t9) + (-t42 * mrSges(5,3) + mrSges(5,1) * t195 + t51 / 0.2e1 + t49 / 0.2e1 - t97 / 0.2e1 + t52 / 0.2e1 + t50 / 0.2e1 - t98 / 0.2e1 + Ifges(5,4) * t99 - t191 * t56 + t192 * t55 - (Ifges(5,2) + t190) * t100 + t147) * t113 + (-mrSges(5,2) * t195 - t17 * t164 + Ifges(5,1) * t99 + (-mrSges(5,3) - t166) * t43 + t174 * mrSges(7,3) + t171 * mrSges(6,3) + (mrSges(6,3) * t150 + mrSges(7,3) * t170 + t165 * t28 + t167 * t63 + t233 * t260 + t237 * t258 + t240 * t256 + t242 * t257) * qJD(5) + (-t157 / 0.2e1 - t159 / 0.2e1) * t56 + (-t163 / 0.2e1 - t161 / 0.2e1) * t55 - t252 * t141 / 0.2e1 + (qJD(5) * t251 + t253) * t235 - (-t153 / 0.2e1 - t155 / 0.2e1 + Ifges(5,4)) * t100) * t112 + (t100 * t88 - t87 * t99) * mrSges(5,3) + t128 * t182 + t19 * t37 + t29 * t38 + t24 * t39 + t30 * t40 + t25 * t45 + t61 * t20 + t4 * t66 + t7 * t67 + t3 * t68 + t8 * t69 + t87 * t21 + t57 * t94 + 0.2e1 * qJD(3) * t250 * qJD(1) + (-t122 * t198 - t139 * t178) * t216; (-mrSges(3,3) - t230) * t142 + (-mrSges(5,3) * t99 + t20 + t21) * t112 + t194 * t109 + (-t140 * t217 + t141 * t218 + t94) * t110 + m(5) * (t109 * t71 + t110 * t72 + t208) + m(7) * (-t109 * t28 + t11 * t200 + t112 * t17 - t201 * t9) + m(6) * (-t109 * t63 + t200 * t23 - t201 * t22 + t208) + (t100 * mrSges(5,3) + t220 * t141 + t221 * t140 + t148 * qJD(5) + m(5) * t42 + m(7) * (-t11 * t197 - t196 * t9 + t173) + m(6) * (-t196 * t22 - t197 * t23 + t172)) * t113 + (-m(4) * t132 - m(5) * t117 - m(7) * t169 - t198 * t216 + t148 - t206 - t254) * qJD(1); t107 * t94 + (m(4) + m(5)) * t195 - t142 * t250 - t194 * t108 + (t218 * t255 - t221) * t141 + (-t217 * t255 + t220) * t140 - m(5) * (-t107 * t72 + t108 * t71) + m(4) * t122 * t178 + t182 + (t108 * t28 + t170 * t255 + t174) * m(7) + (t108 * t63 + t150 * t255 + t171) * m(6); -t219 * t72 - m(7) * (t11 * t18 + t12 * t9 + t28 * t44) + (-t18 + t105) * t66 - (t258 / 0.2e1 - Ifges(5,6)) * t100 + (m(6) * t172 - t140 * t38 + t141 * t40 + (-t140 * t67 - t141 * t69 - t254) * qJD(5)) * pkin(8) + t252 * t235 + t253 * t233 + t256 * t55 / 0.2e1 + t257 * t56 / 0.2e1 + t246 * t108 + (-t12 + t106) * t68 + (-t248 + (m(7) * t28 + t45) * t140 * pkin(5)) * qJD(5) + (-pkin(4) * t43 - t22 * t26 - t23 * t27 - t63 * t72) * m(6) + t172 * mrSges(6,3) + t173 * mrSges(7,3) + t17 * t165 + (-mrSges(5,1) + t167) * t43 - t42 * mrSges(5,2) - t44 * t45 - t27 * t67 - t26 * t69 - t71 * t94 - Ifges(5,5) * t99 + t118 * t37 - t119 * t39 + m(7) * (t1 * t118 + t105 * t11 + t106 * t9 - t119 * t2 + t131 * t17) - pkin(4) * t21 + t131 * t20 - (-(t271 + t272) * t108 - t265) * t107; (-t45 * t92 + t37) * pkin(5) + (t11 * t92 + t9 * t91) * mrSges(7,3) + (t22 * t91 + t23 * t92) * mrSges(6,3) + (-(t10 - t9) * t11 + (-t28 * t92 + t1) * pkin(5)) * m(7) - t97 - t98 + t51 + t52 + t50 + t49 + t147 - t10 * t66 - t22 * t67 + t11 * t68 + t23 * t69 - t28 * (mrSges(7,1) * t92 + mrSges(7,2) * t91) - t63 * (mrSges(6,1) * t92 + mrSges(6,2) * t91) + (t270 * t91 - t266) * t241 + t260 * t240 + (-t261 * t92 + t262 * t91) * t238 + (-t268 * t92 + t251 + t267) * t243; -t91 * t66 + t92 * t68 + 0.2e1 * (t17 / 0.2e1 + t11 * t243 + t9 * t240) * m(7) + t20;];
tauc  = t13(:);
