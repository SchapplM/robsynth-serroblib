% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 17:01
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRRP8_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP8_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:58:14
% EndTime: 2019-12-05 16:58:27
% DurationCPUTime: 4.52s
% Computational Cost: add. (2155->386), mult. (5739->521), div. (0->0), fcn. (3634->8), ass. (0->188)
t116 = cos(qJ(3));
t170 = qJD(2) * t116;
t108 = Ifges(4,4) * t170;
t239 = -t108 / 0.2e1;
t237 = Ifges(5,1) + Ifges(6,1);
t228 = Ifges(6,4) + Ifges(5,5);
t238 = qJD(3) / 0.2e1;
t236 = Ifges(5,6) - Ifges(6,6);
t150 = Ifges(4,5) * t238;
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t113 = sin(qJ(3));
t111 = cos(pkin(5));
t173 = qJD(1) * t111;
t153 = t113 * t173;
t114 = sin(qJ(2));
t110 = sin(pkin(5));
t174 = qJD(1) * t110;
t155 = t114 * t174;
t94 = qJD(2) * pkin(7) + t155;
t71 = t116 * t94 + t153;
t53 = qJD(3) * pkin(8) + t71;
t100 = -pkin(3) * t116 - pkin(8) * t113 - pkin(2);
t117 = cos(qJ(2));
t154 = t117 * t174;
t72 = qJD(2) * t100 - t154;
t15 = -t112 * t53 + t115 * t72;
t16 = t112 * t72 + t115 * t53;
t128 = t112 * t16 + t115 * t15;
t184 = Ifges(6,5) * t115;
t132 = Ifges(6,3) * t112 + t184;
t186 = Ifges(5,4) * t115;
t136 = -Ifges(5,2) * t112 + t186;
t141 = mrSges(6,1) * t112 - mrSges(6,3) * t115;
t143 = mrSges(5,1) * t112 + mrSges(5,2) * t115;
t105 = qJD(4) - t170;
t225 = qJD(5) - t15;
t8 = -pkin(4) * t105 + t225;
t9 = qJ(5) * t105 + t16;
t145 = t112 * t9 - t115 * t8;
t70 = -t113 * t94 + t116 * t173;
t52 = -qJD(3) * pkin(3) - t70;
t164 = t115 * qJD(3);
t172 = qJD(2) * t113;
t89 = t112 * t172 - t164;
t90 = qJD(3) * t112 + t115 * t172;
t17 = pkin(4) * t89 - qJ(5) * t90 + t52;
t204 = t115 / 0.2e1;
t206 = t112 / 0.2e1;
t207 = -t112 / 0.2e1;
t212 = t90 / 0.2e1;
t214 = t89 / 0.2e1;
t215 = -t89 / 0.2e1;
t200 = Ifges(6,5) * t89;
t87 = Ifges(5,4) * t89;
t226 = t228 * t105 + t237 * t90 + t200 - t87;
t229 = t105 / 0.2e1;
t185 = Ifges(6,5) * t112;
t187 = Ifges(5,4) * t112;
t231 = t237 * t115 + t185 - t187;
t86 = Ifges(6,5) * t90;
t30 = Ifges(6,6) * t105 + Ifges(6,3) * t89 + t86;
t201 = Ifges(5,4) * t90;
t33 = -Ifges(5,2) * t89 + Ifges(5,6) * t105 + t201;
t220 = t145 * mrSges(6,2) + t128 * mrSges(5,3) - t132 * t214 - t136 * t215 - t17 * t141 - t52 * t143 - t30 * t206 - t33 * t207 - t231 * t212 - (-t112 * t236 + t115 * t228) * t229 - t226 * t204;
t234 = t172 / 0.2e1;
t183 = qJD(2) * pkin(2);
t95 = -t154 - t183;
t235 = -t95 * mrSges(4,2) + t70 * mrSges(4,3) - Ifges(4,1) * t234 - t150 + t220 + t239;
t163 = qJD(2) * qJD(3);
t147 = t113 * t163;
t162 = qJD(3) * qJD(4);
t167 = qJD(4) * t112;
t61 = t115 * t162 + (-t113 * t167 + t116 * t164) * qJD(2);
t166 = qJD(4) * t115;
t168 = qJD(3) * t116;
t62 = t112 * t162 + (t112 * t168 + t113 * t166) * qJD(2);
t233 = (-Ifges(5,4) + Ifges(6,5)) * t62 + t237 * t61 + t228 * t147;
t232 = t237 * t112 - t184 + t186;
t230 = -m(4) * t70 + m(5) * t52;
t149 = -Ifges(4,6) * qJD(3) / 0.2e1;
t224 = t228 * t112 + t236 * t115;
t177 = t110 * t117;
t151 = qJD(2) * t177;
t123 = qJD(1) * (qJD(3) * t111 + t151);
t169 = qJD(3) * t113;
t36 = t116 * t123 - t169 * t94;
t146 = pkin(3) * t113 - pkin(8) * t116;
t93 = t146 * qJD(3);
t65 = (t93 + t155) * qJD(2);
t3 = t112 * t65 + t115 * t36 + t72 * t166 - t167 * t53;
t4 = -qJD(4) * t16 - t112 * t36 + t115 * t65;
t223 = -t112 * t4 + t115 * t3;
t1 = qJ(5) * t147 + qJD(5) * t105 + t3;
t2 = -pkin(4) * t147 - t4;
t222 = t1 * t115 + t112 * t2;
t165 = qJD(4) * t116;
t25 = pkin(7) * (t112 * t169 - t115 * t165) - t100 * t167 + t115 * t93;
t158 = -Ifges(5,3) / 0.2e1 - Ifges(6,2) / 0.2e1;
t159 = Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1;
t160 = Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1;
t219 = t158 * t105 - t159 * t89 - t160 * t90 - t15 * mrSges(5,1) - t9 * mrSges(6,3) - t95 * mrSges(4,1) - Ifges(5,6) * t215 - Ifges(6,6) * t214 - t149 + (Ifges(4,4) * t113 + t116 * Ifges(4,2)) * qJD(2) / 0.2e1 + t16 * mrSges(5,2) + t71 * mrSges(4,3) + t8 * mrSges(6,1) - (Ifges(5,3) + Ifges(6,2)) * t229 - t228 * t212;
t218 = t61 / 0.2e1;
t217 = -t62 / 0.2e1;
t216 = t62 / 0.2e1;
t213 = -t90 / 0.2e1;
t209 = -t105 / 0.2e1;
t205 = -t115 / 0.2e1;
t203 = mrSges(5,3) * t89;
t202 = mrSges(5,3) * t90;
t37 = t113 * t123 + t168 * t94;
t178 = t110 * t114;
t77 = -t111 * t116 + t113 * t178;
t195 = t37 * t77;
t39 = mrSges(5,1) * t147 - mrSges(5,3) * t61;
t40 = -mrSges(6,1) * t147 + t61 * mrSges(6,2);
t191 = -t39 + t40;
t38 = -mrSges(6,2) * t62 + mrSges(6,3) * t147;
t41 = -mrSges(5,2) * t147 - mrSges(5,3) * t62;
t190 = t41 + t38;
t92 = t146 * qJD(2);
t27 = t112 * t92 + t115 * t70;
t66 = -mrSges(5,2) * t105 - t203;
t69 = -mrSges(6,2) * t89 + mrSges(6,3) * t105;
t189 = t66 + t69;
t67 = mrSges(5,1) * t105 - t202;
t68 = -mrSges(6,1) * t105 + mrSges(6,2) * t90;
t188 = t67 - t68;
t74 = t115 * t116 * pkin(7) + t112 * t100;
t182 = -qJD(3) * mrSges(4,1) + mrSges(5,1) * t89 + mrSges(5,2) * t90 + mrSges(4,3) * t172;
t179 = t100 * t115;
t176 = t112 * t116;
t175 = t115 * t117;
t171 = qJD(2) * t114;
t43 = mrSges(6,1) * t89 - mrSges(6,3) * t90;
t161 = t43 + t182;
t157 = mrSges(4,3) * t170;
t156 = t112 * t177;
t152 = t110 * t171;
t102 = -qJD(3) * mrSges(4,2) + t157;
t148 = -m(4) * t71 - t102;
t144 = mrSges(5,1) * t115 - mrSges(5,2) * t112;
t142 = mrSges(6,1) * t115 + mrSges(6,3) * t112;
t135 = Ifges(5,2) * t115 + t187;
t131 = -Ifges(6,3) * t115 + t185;
t130 = pkin(4) * t115 + qJ(5) * t112;
t129 = pkin(4) * t112 - qJ(5) * t115;
t26 = -t112 * t70 + t115 * t92;
t126 = pkin(7) + t129;
t78 = t111 * t113 + t116 * t178;
t47 = t110 * t175 + t112 * t78;
t122 = -t4 * mrSges(5,1) + t2 * mrSges(6,1) + t3 * mrSges(5,2) - t1 * mrSges(6,3);
t24 = t112 * t93 + t100 * t166 + (-t112 * t165 - t113 * t164) * pkin(7);
t118 = qJD(2) ^ 2;
t104 = Ifges(6,2) * t147;
t103 = Ifges(5,3) * t147;
t96 = -pkin(3) - t130;
t91 = (-mrSges(4,1) * t116 + mrSges(4,2) * t113) * qJD(2);
t83 = (mrSges(4,1) * t113 + mrSges(4,2) * t116) * t163;
t76 = qJD(4) * t129 - qJD(5) * t112;
t75 = t126 * t113;
t73 = -pkin(7) * t176 + t179;
t64 = (t112 * t114 + t116 * t175) * t174;
t63 = -t115 * t155 + t154 * t176;
t60 = -t179 + (pkin(7) * t112 + pkin(4)) * t116;
t59 = -qJ(5) * t116 + t74;
t58 = Ifges(6,4) * t61;
t57 = Ifges(5,5) * t61;
t56 = Ifges(5,6) * t62;
t55 = Ifges(6,6) * t62;
t48 = t115 * t78 - t156;
t46 = -qJD(3) * t77 + t116 * t151;
t45 = qJD(3) * t78 + t113 * t151;
t42 = pkin(4) * t90 + qJ(5) * t89;
t29 = t153 + (qJD(2) * t129 + t94) * t116;
t23 = (qJD(4) * t130 - qJD(5) * t115) * t113 + t126 * t168;
t22 = -pkin(4) * t172 - t26;
t21 = qJ(5) * t172 + t27;
t20 = -pkin(4) * t169 - t25;
t19 = mrSges(5,1) * t62 + mrSges(5,2) * t61;
t18 = mrSges(6,1) * t62 - mrSges(6,3) * t61;
t14 = qJ(5) * t169 - qJD(5) * t116 + t24;
t11 = t61 * Ifges(5,4) - t62 * Ifges(5,2) + Ifges(5,6) * t147;
t10 = t61 * Ifges(6,5) + Ifges(6,6) * t147 + t62 * Ifges(6,3);
t7 = -qJD(4) * t47 + t112 * t152 + t115 * t46;
t6 = -qJD(4) * t156 + t112 * t46 - t115 * t152 + t166 * t78;
t5 = pkin(4) * t62 - qJ(5) * t61 - qJD(5) * t90 + t37;
t12 = [-t78 * mrSges(4,3) * t147 + t46 * t102 + t189 * t7 - t188 * t6 + t190 * t48 + t191 * t47 + (qJD(3) * t157 + t18 + t19) * t77 + t161 * t45 + ((-mrSges(3,2) * t118 - t83) * t117 + (-mrSges(3,1) * t118 + qJD(2) * t91) * t114) * t110 + m(4) * (t36 * t78 + t195 - t45 * t70 + t46 * t71 + (t95 - t154) * t152) + m(5) * (-t15 * t6 + t16 * t7 + t3 * t48 - t4 * t47 + t45 * t52 + t195) + m(6) * (t1 * t48 + t17 * t45 + t2 * t47 + t5 * t77 + t6 * t8 + t7 * t9); -t91 * t155 - pkin(2) * t83 + t14 * t69 + t75 * t18 + t20 * t68 + t23 * t43 + t24 * t66 + t25 * t67 + t59 * t38 + t73 * t39 + t60 * t40 + t74 * t41 - t189 * t64 + t188 * t63 - m(6) * (t63 * t8 + t64 * t9) - m(5) * (-t15 * t63 + t16 * t64) + m(5) * (t15 * t25 + t16 * t24 + t3 * t74 + t4 * t73) + m(6) * (t1 * t59 + t14 * t9 + t17 * t23 + t2 * t60 + t20 * t8 + t5 * t75) + 0.2e1 * (-t95 / 0.2e1 - t183 / 0.2e1) * m(4) * t155 + (-t104 / 0.2e1 - t103 / 0.2e1 - t58 / 0.2e1 - t55 / 0.2e1 + t56 / 0.2e1 - t57 / 0.2e1 - t159 * t62 - t160 * t61 + (m(4) * pkin(7) + mrSges(4,3)) * t36 + (-mrSges(4,1) * t171 + t117 * t148) * t174 + ((t182 + t230) * pkin(7) + 0.3e1 / 0.2e1 * t108 + t150 - t235) * qJD(3) + t122) * t116 + (t11 * t207 + t132 * t216 + t136 * t217 + t5 * t141 + t10 * t206 + (-t112 * t3 - t115 * t4) * mrSges(5,3) + (-t1 * t112 + t115 * t2) * mrSges(6,2) + (t149 - t219) * qJD(3) + (t131 * t215 + t135 * t214 + t52 * t144 + t17 * t142 + t33 * t205 + (t112 * t15 - t115 * t16) * mrSges(5,3) + (-t112 * t8 - t115 * t9) * mrSges(6,2) + t232 * t213 + t224 * t209 + t226 * t207) * qJD(4) + (mrSges(4,2) * t155 + ((-0.3e1 / 0.2e1 * Ifges(4,4) + t160 * t115 + t159 * t112) * t113 + (0.3e1 / 0.2e1 * Ifges(4,1) - 0.3e1 / 0.2e1 * Ifges(4,2) + t158) * t116) * qJD(3)) * qJD(2) + t231 * t218 + (qJD(4) * t30 + t233) * t204 + (t148 * qJD(3) + t19) * pkin(7) + (-m(6) * t17 - t161 - t230) * t154 + (mrSges(4,3) + t143 + (m(4) + m(5)) * pkin(7)) * t37) * t113; -t182 * t71 + t222 * mrSges(6,2) + m(6) * (pkin(8) * t222 + t17 * t76 + t5 * t96) + t223 * mrSges(5,3) + m(5) * (-pkin(3) * t37 + pkin(8) * t223) + t96 * t18 - t70 * t102 - t22 * t68 - t21 * t69 - t27 * t66 - t26 * t67 - t36 * mrSges(4,2) - pkin(3) * t19 - t5 * t142 + (-mrSges(4,1) - t144) * t37 + (t112 * t191 + t115 * t190) * pkin(8) + (t76 - t29) * t43 - m(5) * (t15 * t26 + t16 * t27 + t52 * t71) - m(6) * (t17 * t29 + t21 * t9 + t22 * t8) + t131 * t216 + t135 * t217 + t11 * t204 + t10 * t205 + t233 * t206 + t232 * t218 + ((-m(5) * t128 - m(6) * t145 - t112 * t189 - t115 * t188) * pkin(8) - t220) * qJD(4) + ((Ifges(4,4) * t234 + t224 * t238 + t149 + t219) * t113 + (t239 + (Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1) * t172 + t150 + t235) * t116) * qJD(2); (t8 * t89 + t9 * t90) * mrSges(6,2) + t104 + t103 + (t188 + t202) * t16 + (-t189 - t203) * t15 + t58 + t55 - t56 + t57 - t52 * (mrSges(5,1) * t90 - mrSges(5,2) * t89) - t17 * (t90 * mrSges(6,1) + t89 * mrSges(6,3)) + qJD(5) * t69 + qJ(5) * t38 - pkin(4) * t40 - t42 * t43 - t122 + t33 * t212 + (Ifges(6,3) * t90 - t200) * t215 + (-t228 * t89 - t236 * t90) * t209 + (-pkin(4) * t2 + qJ(5) * t1 - t16 * t8 - t17 * t42 + t225 * t9) * m(6) + (-Ifges(5,2) * t90 + t226 - t87) * t214 + (-t237 * t89 - t201 + t30 + t86) * t213; -t105 * t69 + t90 * t43 + 0.2e1 * (t2 / 0.2e1 + t9 * t209 + t17 * t212) * m(6) + t40;];
tauc = t12(:);
