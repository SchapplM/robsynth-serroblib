% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRRR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRRR13_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:33
% EndTime: 2019-12-31 19:14:43
% DurationCPUTime: 4.71s
% Computational Cost: add. (3704->406), mult. (8262->585), div. (0->0), fcn. (4954->6), ass. (0->201)
t139 = sin(qJ(4));
t140 = sin(qJ(3));
t137 = t140 * qJD(1);
t178 = t139 * t137;
t238 = -pkin(8) - pkin(7);
t179 = qJD(4) * t238;
t143 = cos(qJ(3));
t169 = pkin(3) * t143 + pkin(7) * t140;
t120 = t169 * qJD(1);
t144 = -pkin(1) - pkin(6);
t132 = qJD(1) * t144 + qJD(2);
t142 = cos(qJ(4));
t193 = t142 * t143;
t68 = t139 * t120 + t132 * t193;
t274 = -pkin(8) * t178 + t139 * t179 - t68;
t182 = pkin(8) * t140 * t142;
t195 = t139 * t143;
t67 = t142 * t120 - t132 * t195;
t273 = t142 * t179 - (pkin(4) * t143 + t182) * qJD(1) - t67;
t196 = t132 * t143;
t107 = -qJD(3) * pkin(3) - t196;
t190 = qJD(3) * t142;
t192 = qJD(1) * t143;
t115 = -t139 * t192 + t190;
t135 = t137 + qJD(4);
t124 = pkin(3) * t140 - pkin(7) * t143 + qJ(2);
t103 = t124 * qJD(1);
t123 = t140 * t132;
t106 = qJD(3) * pkin(7) + t123;
t55 = t142 * t103 - t106 * t139;
t56 = t103 * t139 + t106 * t142;
t159 = t139 * t56 + t142 * t55;
t214 = Ifges(5,4) * t142;
t164 = -Ifges(5,2) * t139 + t214;
t215 = Ifges(5,4) * t139;
t166 = Ifges(5,1) * t142 - t215;
t167 = mrSges(5,1) * t139 + mrSges(5,2) * t142;
t212 = Ifges(5,6) * t139;
t213 = Ifges(5,5) * t142;
t226 = t142 / 0.2e1;
t227 = -t139 / 0.2e1;
t116 = qJD(3) * t139 + t142 * t192;
t232 = t116 / 0.2e1;
t216 = Ifges(5,4) * t116;
t52 = t115 * Ifges(5,2) + t135 * Ifges(5,6) + t216;
t109 = Ifges(5,4) * t115;
t53 = t116 * Ifges(5,1) + t135 * Ifges(5,5) + t109;
t272 = -t159 * mrSges(5,3) + t107 * t167 + (-t212 + t213) * t135 / 0.2e1 + t164 * t115 / 0.2e1 + t166 * t232 + t226 * t53 + t227 * t52;
t138 = sin(qJ(5));
t141 = cos(qJ(5));
t46 = pkin(8) * t115 + t56;
t203 = t141 * t46;
t45 = -pkin(8) * t116 + t55;
t37 = pkin(4) * t135 + t45;
t10 = t138 * t37 + t203;
t189 = qJD(3) * t143;
t172 = qJD(1) * t189;
t133 = Ifges(6,3) * t172;
t170 = t141 * t115 - t116 * t138;
t64 = t115 * t138 + t116 * t141;
t225 = Ifges(6,4) * t64;
t131 = qJD(5) + t135;
t231 = -t131 / 0.2e1;
t245 = -t64 / 0.2e1;
t247 = -t170 / 0.2e1;
t57 = Ifges(6,4) * t170;
t25 = Ifges(6,1) * t64 + Ifges(6,5) * t131 + t57;
t71 = -pkin(4) * t115 + t107;
t205 = t138 * t46;
t9 = t141 * t37 - t205;
t271 = t133 + (Ifges(6,5) * t170 - Ifges(6,6) * t64) * t231 + (t10 * t64 + t170 * t9) * mrSges(6,3) + (-Ifges(6,2) * t64 + t25 + t57) * t247 - t71 * (mrSges(6,1) * t64 + mrSges(6,2) * t170) + (Ifges(6,1) * t170 - t225) * t245;
t128 = t238 * t139;
t129 = t238 * t142;
t74 = t128 * t141 + t129 * t138;
t270 = qJD(5) * t74 + t273 * t138 + t274 * t141;
t75 = t128 * t138 - t129 * t141;
t269 = -qJD(5) * t75 - t274 * t138 + t273 * t141;
t155 = t138 * t139 - t141 * t142;
t255 = qJD(4) + qJD(5);
t268 = t255 * t155;
t200 = Ifges(4,5) * qJD(3);
t218 = Ifges(4,4) * t140;
t267 = t200 / 0.2e1 + (t143 * Ifges(4,1) - t218) * qJD(1) / 0.2e1 + t272;
t183 = qJD(3) * qJD(4);
t185 = qJD(4) * t143;
t76 = t142 * t183 + (-t139 * t185 - t140 * t190) * qJD(1);
t191 = qJD(3) * t140;
t150 = t139 * t191 - t142 * t185;
t77 = qJD(1) * t150 - t139 * t183;
t19 = t170 * qJD(5) + t138 * t77 + t141 * t76;
t113 = qJD(3) * t169 + qJD(2);
t90 = t113 * qJD(1);
t149 = -qJD(4) * t56 + t142 * t90;
t13 = -pkin(8) * t76 + (pkin(4) * qJD(1) - t132 * t139) * t189 + t149;
t176 = t132 * t189;
t186 = qJD(4) * t142;
t187 = qJD(4) * t139;
t26 = t103 * t186 - t106 * t187 + t139 * t90 + t142 * t176;
t14 = pkin(8) * t77 + t26;
t2 = qJD(5) * t9 + t13 * t138 + t14 * t141;
t20 = -qJD(5) * t64 - t138 * t76 + t141 * t77;
t3 = -qJD(5) * t10 + t13 * t141 - t138 * t14;
t266 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20;
t265 = (qJ(2) * (m(3) + m(4)) + mrSges(3,3)) * qJD(1);
t24 = Ifges(6,2) * t170 + Ifges(6,6) * t131 + t225;
t263 = t24 / 0.2e1;
t173 = -Ifges(4,6) * qJD(3) / 0.2e1;
t118 = t138 * t142 + t139 * t141;
t92 = t118 * t143;
t258 = t155 * qJD(1) - qJD(3) * t92 + t140 * t268;
t101 = t118 * qJD(1);
t66 = t255 * t118;
t94 = t155 * t143;
t257 = -qJD(3) * t94 - t140 * t66 - t101;
t31 = -mrSges(6,1) * t170 + mrSges(6,2) * t64;
t180 = m(6) * t71 + t31;
t256 = t155 * t140;
t194 = t140 * t144;
t130 = t142 * t194;
t82 = t139 * t124 + t130;
t27 = -t139 * t176 + t149;
t160 = -t139 * t27 + t142 * t26;
t79 = -mrSges(5,2) * t135 + mrSges(5,3) * t115;
t80 = mrSges(5,1) * t135 - mrSges(5,3) * t116;
t253 = -m(5) * t159 - t139 * t79 - t142 * t80;
t251 = t27 * mrSges(5,1) - t26 * mrSges(5,2) + Ifges(5,5) * t76 + Ifges(5,6) * t77 + t266;
t249 = t19 / 0.2e1;
t248 = t20 / 0.2e1;
t246 = t170 / 0.2e1;
t244 = t64 / 0.2e1;
t243 = t76 / 0.2e1;
t242 = t77 / 0.2e1;
t241 = -t92 / 0.2e1;
t240 = -t94 / 0.2e1;
t235 = -t115 / 0.2e1;
t233 = -t116 / 0.2e1;
t230 = t131 / 0.2e1;
t229 = -t135 / 0.2e1;
t224 = pkin(4) * t139;
t86 = qJD(1) * t256;
t220 = -t268 - t86;
t85 = t140 * t101;
t219 = -t66 - t85;
t217 = Ifges(4,4) * t143;
t211 = qJ(2) * mrSges(4,1);
t210 = qJ(2) * mrSges(4,2);
t201 = -qJD(3) * mrSges(4,1) - mrSges(5,1) * t115 + mrSges(5,2) * t116 + mrSges(4,3) * t192;
t198 = qJD(3) * mrSges(4,2);
t188 = qJD(3) * t144;
t184 = qJD(1) * qJD(2);
t181 = t139 * t194;
t177 = t132 * t191;
t175 = t143 * t188;
t174 = -t200 / 0.2e1;
t171 = -t139 * t144 + pkin(4);
t168 = mrSges(5,1) * t142 - mrSges(5,2) * t139;
t165 = Ifges(5,1) * t139 + t214;
t163 = Ifges(5,2) * t142 + t215;
t161 = Ifges(5,5) * t139 + Ifges(5,6) * t142;
t112 = t142 * t124;
t60 = -pkin(8) * t193 + t140 * t171 + t112;
t70 = -pkin(8) * t195 + t82;
t29 = -t138 * t70 + t141 * t60;
t30 = t138 * t60 + t141 * t70;
t158 = t139 * t55 - t142 * t56;
t58 = mrSges(5,1) * t172 - mrSges(5,3) * t76;
t59 = -mrSges(5,2) * t172 + mrSges(5,3) * t77;
t157 = -t139 * t58 + t142 * t59;
t43 = -qJD(4) * t181 + t139 * t113 + t124 * t186 + t142 * t175;
t147 = t55 * mrSges(5,1) + t9 * mrSges(6,1) + t131 * Ifges(6,3) + t64 * Ifges(6,5) + t170 * Ifges(6,6) + t135 * Ifges(5,3) + t116 * Ifges(5,5) + t115 * Ifges(5,6) + t173 - (-Ifges(4,2) * t140 + t217) * qJD(1) / 0.2e1 - t10 * mrSges(6,2) - t56 * mrSges(5,2);
t136 = -pkin(4) * t142 - pkin(3);
t134 = Ifges(5,3) * t172;
t126 = -mrSges(4,3) * t137 - t198;
t119 = qJD(1) * (t140 * mrSges(4,1) + t143 * mrSges(4,2));
t114 = (-t144 + t224) * t143;
t96 = t142 * t113;
t91 = t118 * t140;
t87 = -pkin(4) * t178 + t123;
t81 = t112 - t181;
t78 = -pkin(4) * t150 + t140 * t188;
t50 = -pkin(4) * t77 + t177;
t48 = mrSges(6,1) * t131 - mrSges(6,3) * t64;
t47 = -mrSges(6,2) * t131 + mrSges(6,3) * t170;
t44 = -t82 * qJD(4) - t139 * t175 + t96;
t42 = -mrSges(5,1) * t77 + mrSges(5,2) * t76;
t39 = t76 * Ifges(5,1) + t77 * Ifges(5,4) + Ifges(5,5) * t172;
t38 = t76 * Ifges(5,4) + t77 * Ifges(5,2) + Ifges(5,6) * t172;
t36 = t118 * t191 + t143 * t268;
t34 = qJD(3) * t256 - t143 * t66;
t32 = pkin(8) * t150 + t43;
t28 = t96 + (-t130 + (pkin(8) * t143 - t124) * t139) * qJD(4) + (t143 * t171 + t182) * qJD(3);
t16 = -mrSges(6,2) * t172 + mrSges(6,3) * t20;
t15 = mrSges(6,1) * t172 - mrSges(6,3) * t19;
t12 = t141 * t45 - t205;
t11 = -t138 * t45 - t203;
t8 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t7 = t19 * Ifges(6,1) + t20 * Ifges(6,4) + Ifges(6,5) * t172;
t6 = t19 * Ifges(6,4) + t20 * Ifges(6,2) + Ifges(6,6) * t172;
t5 = -qJD(5) * t30 - t138 * t32 + t141 * t28;
t4 = qJD(5) * t29 + t138 * t28 + t141 * t32;
t1 = [m(6) * (t10 * t4 + t114 * t50 + t2 * t30 + t29 * t3 + t5 * t9 + t71 * t78) + t114 * t8 + t78 * t31 + t43 * t79 + t44 * t80 + t81 * t58 + t82 * t59 + t71 * (-mrSges(6,1) * t36 + mrSges(6,2) * t34) + m(5) * (t26 * t82 + t27 * t81 + t43 * t56 + t44 * t55) + t4 * t47 + t5 * t48 + t29 * t15 + t30 * t16 + t34 * t25 / 0.2e1 + t50 * (mrSges(6,1) * t92 - mrSges(6,2) * t94) + (t10 * t36 - t2 * t92 + t3 * t94 - t34 * t9) * mrSges(6,3) + (-Ifges(6,4) * t94 - Ifges(6,2) * t92) * t248 + (-Ifges(6,1) * t94 - Ifges(6,4) * t92) * t249 + (t119 + 0.2e1 * t265) * qJD(2) + (t133 / 0.2e1 + t134 / 0.2e1 + mrSges(4,1) * t184 + (t174 + (-0.2e1 * t210 + 0.3e1 / 0.2e1 * t218) * qJD(1) + (m(5) * t107 + t201) * t144 - t267) * qJD(3) + t251) * t140 + (Ifges(6,4) * t34 + Ifges(6,2) * t36) * t246 + (Ifges(6,5) * t34 + Ifges(6,6) * t36) * t230 + t7 * t240 + t6 * t241 + (Ifges(6,1) * t34 + Ifges(6,4) * t36) * t244 + t36 * t263 + (t39 * t226 + t38 * t227 - t144 * t42 + t166 * t243 + t164 * t242 + mrSges(4,2) * t184 + (-t139 * t26 - t142 * t27) * mrSges(5,3) + (t161 * t229 + t163 * t235 + t165 * t233 + t107 * t168 - t142 * t52 / 0.2e1 + t53 * t227 + t158 * mrSges(5,3)) * qJD(4) + (t147 + (-m(5) * t144 + t167) * t123 + (0.2e1 * t211 + (t213 / 0.2e1 - t212 / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,4)) * t143 + Ifges(6,5) * t240 + Ifges(6,6) * t241 + (Ifges(5,3) / 0.2e1 + 0.3e1 / 0.2e1 * Ifges(4,2) + Ifges(6,3) / 0.2e1 - 0.3e1 / 0.2e1 * Ifges(4,1)) * t140) * qJD(1) + t144 * t126 + t173) * qJD(3)) * t143; -t91 * t15 - t256 * t16 + t258 * t48 + t257 * t47 + (-t42 - t8 + (-t139 * t80 + t142 * t79 + t126) * qJD(3)) * t143 - m(5) * t158 * t189 + (-t80 * t186 - t79 * t187 + m(5) * (-t55 * t186 - t56 * t187 + t160) + (m(5) * (t107 - t196) + t180 + t201) * qJD(3) + t157) * t140 + (t257 * t10 - t143 * t50 - t2 * t256 + t258 * t9 - t3 * t91) * m(6) + (-t119 + t253 - t265) * qJD(1); t157 * pkin(7) + t160 * mrSges(5,3) + (t253 * pkin(7) + t180 * t224 + t272) * qJD(4) + (-pkin(3) * t177 + t160 * pkin(7) - t107 * t123 - t55 * t67 - t56 * t68) * m(5) + t139 * t39 / 0.2e1 + t136 * t8 + t118 * t7 / 0.2e1 - t87 * t31 - t68 * t79 - t67 * t80 + t74 * t15 + t75 * t16 - pkin(3) * t42 + t269 * t48 + t270 * t47 + (t270 * t10 + t136 * t50 + t2 * t75 + t269 * t9 + t3 * t74 - t71 * t87) * m(6) + ((t174 + (t210 - t218 / 0.2e1) * qJD(1) + t267) * t140 + (-t147 + (-t211 + t217 / 0.2e1 + (Ifges(4,1) / 0.2e1 - Ifges(4,2) / 0.2e1) * t140) * qJD(1) + t173 + (Ifges(6,5) * t118 - Ifges(6,6) * t155 + t161) * qJD(3) / 0.2e1) * t143) * qJD(1) - t155 * t6 / 0.2e1 + t50 * (mrSges(6,1) * t155 + mrSges(6,2) * t118) + (Ifges(6,4) * t118 - Ifges(6,2) * t155) * t248 + (Ifges(6,1) * t118 - Ifges(6,4) * t155) * t249 + (t10 * t219 - t118 * t3 - t155 * t2 - t220 * t9) * mrSges(6,3) + (-t85 / 0.2e1 - t66 / 0.2e1) * t24 + (Ifges(6,4) * t86 + Ifges(6,2) * t85) * t247 + (Ifges(6,5) * t86 + Ifges(6,6) * t85) * t231 + t163 * t242 + t165 * t243 + (Ifges(6,1) * t86 + Ifges(6,4) * t85) * t245 + t38 * t226 + ((-t126 - t198) * t143 + ((-mrSges(4,1) - t168) * qJD(3) - t201) * t140) * t132 + (-mrSges(6,1) * t219 + mrSges(6,2) * t220) * t71 + (-t86 / 0.2e1 - t268 / 0.2e1) * t25 + (-Ifges(6,4) * t268 - Ifges(6,2) * t66) * t246 + (-Ifges(6,5) * t268 - Ifges(6,6) * t66) * t230 + (-Ifges(6,1) * t268 - Ifges(6,4) * t66) * t244; (t138 * t16 + t141 * t15 + m(6) * (t138 * t2 + t141 * t3) - t180 * t116 + (-t138 * t48 + t141 * t47 + m(6) * (t10 * t141 - t138 * t9)) * qJD(5)) * pkin(4) + t251 + (-Ifges(5,2) * t116 + t109 + t53) * t235 - m(6) * (t10 * t12 + t11 * t9) + t134 - t107 * (mrSges(5,1) * t116 + mrSges(5,2) * t115) - t55 * t79 + t56 * t80 - t12 * t47 - t11 * t48 + t64 * t263 + (Ifges(5,5) * t115 - Ifges(5,6) * t116) * t229 + t52 * t232 + (Ifges(5,1) * t115 - t216) * t233 + (t115 * t55 + t116 * t56) * mrSges(5,3) + t271; t10 * t48 + t24 * t244 - t9 * t47 + t266 + t271;];
tauc = t1(:);
