% Calculate vector of inverse dynamics joint torques for
% S5RRPRP1
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,theta3]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:22
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP1_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP1_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:22:00
% EndTime: 2019-12-05 18:22:07
% DurationCPUTime: 2.72s
% Computational Cost: add. (2026->325), mult. (3250->405), div. (0->0), fcn. (1649->12), ass. (0->151)
t221 = Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1;
t220 = mrSges(5,2) + mrSges(6,2);
t219 = -mrSges(5,3) - mrSges(6,3);
t218 = Ifges(5,1) + Ifges(6,1);
t216 = Ifges(6,5) + Ifges(5,5);
t215 = Ifges(5,2) + Ifges(6,2);
t214 = Ifges(6,6) + Ifges(5,6);
t135 = -qJ(5) - pkin(7);
t213 = m(6) * t135 + mrSges(4,2);
t130 = qJDD(1) + qJDD(2);
t133 = sin(pkin(8));
t134 = cos(pkin(8));
t137 = sin(qJ(2));
t189 = qJD(1) * pkin(1);
t169 = t137 * t189;
t140 = cos(qJ(2));
t202 = pkin(1) * t140;
t85 = -qJD(2) * t169 + qJDD(1) * t202;
t63 = pkin(2) * t130 + t85;
t176 = qJD(2) * t140;
t86 = (qJD(1) * t176 + qJDD(1) * t137) * pkin(1);
t26 = t133 * t63 + t134 * t86;
t18 = pkin(7) * t130 + t26;
t212 = qJD(3) * qJD(4) + t18;
t139 = cos(qJ(4));
t211 = (Ifges(5,4) + Ifges(6,4)) * t139;
t136 = sin(qJ(4));
t210 = t136 * t220;
t197 = mrSges(5,1) + mrSges(6,1);
t175 = qJD(3) * t136;
t131 = qJD(1) + qJD(2);
t168 = t140 * t189;
t93 = pkin(2) * t131 + t168;
t45 = t133 * t93 + t134 * t169;
t33 = pkin(7) * t131 + t45;
t29 = t139 * t33 + t175;
t178 = t29 * qJD(4);
t123 = t139 * qJDD(3);
t6 = -t136 * t18 + t123 - t178;
t208 = -t6 - t178;
t201 = pkin(2) * t133;
t200 = pkin(2) * t134;
t199 = pkin(4) * t139;
t174 = qJD(4) * t136;
t5 = t136 * qJDD(3) + t139 * t212 - t174 * t33;
t198 = t139 * t5;
t182 = t131 * t136;
t89 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t182;
t90 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t182;
t195 = t89 + t90;
t181 = t131 * t139;
t91 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t181;
t92 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t181;
t194 = t91 + t92;
t193 = Ifges(5,4) * t136;
t191 = Ifges(6,4) * t136;
t118 = pkin(2) + t202;
t179 = t134 * t137;
t72 = pkin(1) * t179 + t133 * t118;
t65 = pkin(7) + t72;
t188 = t139 * t65;
t180 = t133 * t137;
t69 = (t134 * t140 - t180) * qJD(2) * pkin(1);
t187 = t139 * t69;
t185 = -qJ(5) - t65;
t132 = qJ(1) + qJ(2);
t114 = pkin(7) + t201;
t177 = -qJ(5) - t114;
t173 = qJD(4) * t139;
t125 = t139 * qJD(5);
t172 = m(4) + m(5) + m(6);
t167 = pkin(4) * t174;
t117 = pkin(3) + t199;
t166 = t131 * t174;
t77 = t130 * t139 - t166;
t78 = t130 * t136 + t131 * t173;
t30 = -t77 * mrSges(6,1) + t78 * mrSges(6,2);
t25 = -t133 * t86 + t134 * t63;
t163 = qJ(5) * t131 + t33;
t102 = t133 * t169;
t44 = t134 * t93 - t102;
t162 = qJD(4) * t185;
t71 = -pkin(1) * t180 + t118 * t134;
t161 = qJD(4) * t177;
t64 = -pkin(3) - t71;
t126 = t139 * qJD(3);
t15 = -t136 * t163 + t126;
t13 = qJD(4) * pkin(4) + t15;
t28 = -t136 * t33 + t126;
t157 = t28 * mrSges(5,3) + t13 * mrSges(6,3);
t16 = t139 * t163 + t175;
t156 = t29 * mrSges(5,3) + t16 * mrSges(6,3);
t124 = pkin(8) + t132;
t112 = sin(t124);
t113 = cos(t124);
t155 = g(2) * t113 + g(3) * t112;
t154 = pkin(2) * t172 + mrSges(3,1);
t100 = -t139 * mrSges(5,1) + mrSges(5,2) * t136;
t153 = -t139 * mrSges(6,1) + mrSges(6,2) * t136;
t152 = Ifges(5,2) * t139 + t193;
t151 = Ifges(6,2) * t139 + t191;
t150 = -t13 * t136 + t139 * t16;
t149 = -t136 * t28 + t139 * t29;
t148 = m(5) * pkin(3) + m(6) * t117 + mrSges(4,1);
t147 = mrSges(2,1) + (m(3) + t172) * pkin(1);
t17 = -pkin(3) * t130 - t25;
t146 = pkin(1) * (t133 * t140 + t179);
t67 = qJD(2) * t146;
t145 = t139 * t197 + t148;
t127 = sin(t132);
t128 = cos(t132);
t144 = t128 * mrSges(3,2) - t112 * t210 + t113 * t219 + t127 * t154;
t143 = -t127 * mrSges(3,2) + t154 * t128 - t113 * t210 + (m(5) * pkin(7) - t213 - t219) * t112;
t27 = -t117 * t131 + qJD(5) - t44;
t3 = qJ(5) * t77 + t125 * t131 + t5;
t32 = -pkin(3) * t131 - t44;
t59 = Ifges(6,6) * qJD(4) + t131 * t151;
t60 = Ifges(5,6) * qJD(4) + t131 * t152;
t103 = Ifges(6,4) * t181;
t61 = Ifges(6,1) * t182 + Ifges(6,5) * qJD(4) + t103;
t104 = Ifges(5,4) * t181;
t62 = Ifges(5,1) * t182 + Ifges(5,5) * qJD(4) + t104;
t8 = -pkin(4) * t77 + qJDD(5) + t17;
t142 = t3 * t139 * mrSges(6,3) + t85 * mrSges(3,1) + t25 * mrSges(4,1) - t86 * mrSges(3,2) + mrSges(5,3) * t198 + t17 * t100 + t8 * t153 + (-t214 * t136 + t216 * t139) * qJD(4) ^ 2 / 0.2e1 - (t60 + t59) * t174 / 0.2e1 + (t139 * t218 - t191 - t193) * t166 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t130 + (t32 * (mrSges(5,1) * t136 + mrSges(5,2) * t139) + t27 * (mrSges(6,1) * t136 + mrSges(6,2) * t139)) * qJD(4) + (t152 / 0.2e1 + t151 / 0.2e1 + t136 * t221 + t215 * t139 / 0.2e1) * t77 + (t218 * t136 + t211 / 0.2e1 + t139 * t221) * t78 + (t62 + t61 + (-t136 * t215 + t211) * t131) * t173 / 0.2e1 + (t136 * t216 + t139 * t214) * qJDD(4);
t141 = cos(qJ(1));
t138 = sin(qJ(1));
t129 = t139 * qJ(5);
t115 = -pkin(3) - t200;
t108 = t113 * pkin(7);
t98 = -t117 - t200;
t84 = t114 * t139 + t129;
t83 = t177 * t136;
t76 = t100 * t131;
t75 = t153 * t131;
t68 = t134 * t168 - t102;
t66 = qJD(1) * t146;
t54 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t78;
t53 = qJDD(4) * mrSges(6,1) - mrSges(6,3) * t78;
t52 = -qJDD(4) * mrSges(5,2) + mrSges(5,3) * t77;
t51 = -qJDD(4) * mrSges(6,2) + mrSges(6,3) * t77;
t50 = -qJD(5) * t136 + t139 * t161;
t49 = t136 * t161 + t125;
t48 = t64 - t199;
t46 = t67 + t167;
t39 = t129 + t188;
t38 = t185 * t136;
t31 = -mrSges(5,1) * t77 + mrSges(5,2) * t78;
t11 = (-qJD(5) - t69) * t136 + t139 * t162;
t10 = t136 * t162 + t125 + t187;
t2 = -t33 * t173 + qJDD(4) * pkin(4) - qJ(5) * t78 + t123 + (-qJD(5) * t131 - t212) * t136;
t1 = [m(5) * (-t28 * t65 * t173 + t17 * t64 + t29 * t187 + t5 * t188 + t32 * t67) + ((-t16 * qJD(4) - t2) * mrSges(6,3) + t208 * mrSges(5,3) + (-m(5) * t28 - t90) * t69 + (m(5) * t208 - qJD(4) * t92 - t54) * t65) * t136 + (-t138 * mrSges(2,2) + t113 * t145 + t141 * t147 + t143) * g(2) + t11 * t89 + t10 * t91 + t46 * t75 + t67 * t76 + t142 + t48 * t30 + t39 * t51 + t38 * t53 + t64 * t31 + m(6) * (t10 * t16 + t11 * t13 + t2 * t38 + t27 * t46 + t3 * t39 + t48 * t8) + m(4) * (t25 * t71 + t26 * t72 - t44 * t67 + t45 * t69) + Ifges(2,3) * qJDD(1) + (t71 * t130 - t67 * t131) * mrSges(4,1) + (-t130 * t72 - t131 * t69 - t26) * mrSges(4,2) + (t65 * t52 + t69 * t92 + (-t65 * t90 - t157) * qJD(4)) * t139 + (-m(5) * t108 + mrSges(2,2) * t141 + t145 * t112 + t113 * t213 + t147 * t138 + t144) * g(3) + (m(3) * (t137 * t86 + t140 * t85) + (-t130 * t137 - t131 * t176) * mrSges(3,2) + (-qJD(2) * t131 * t137 + t130 * t140) * mrSges(3,1)) * pkin(1); t115 * t31 + t98 * t30 + t83 * t53 + t84 * t51 + t50 * t89 + t49 * t91 + t142 + (-t6 * mrSges(5,3) - t2 * mrSges(6,3) - t114 * t54 + t195 * t68 + (pkin(4) * t75 - t114 * t92 - t156) * qJD(4)) * t136 + (t114 * t52 - t194 * t68 + (-t114 * t90 - t157) * qJD(4) + t197 * t155) * t139 + t130 * mrSges(4,1) * t200 + (t113 * t148 + t143) * g(2) + (-t130 * t201 - t26) * mrSges(4,2) + (-t75 - t76) * t66 + (t66 * mrSges(4,1) + t68 * mrSges(4,2) + (mrSges(3,1) * t137 + mrSges(3,2) * t140) * t189) * t131 + (t112 * mrSges(4,1) + t113 * mrSges(4,2) + t144) * g(3) + (-t150 * t68 + t13 * t50 + t16 * t49 + t2 * t83 + t3 * t84 + t8 * t98 + (t112 * t117 + t113 * t135) * g(3) + (-t66 + t167) * t27) * m(6) + (t44 * t66 - t45 * t68 + (t133 * t26 + t134 * t25) * pkin(2)) * m(4) + (-t149 * t68 - t32 * t66 + (-t136 * t6 + t198 + (-t136 * t29 - t139 * t28) * qJD(4)) * t114 + t115 * t17 + (pkin(3) * t112 - t108) * g(3)) * m(5); m(4) * qJDD(3) + (t53 + t54) * t139 + (t51 + t52) * t136 + (-t195 * t136 + t194 * t139) * qJD(4) + m(5) * (qJD(4) * t149 + t136 * t5 + t139 * t6) + m(6) * (qJD(4) * t150 + t136 * t3 + t139 * t2) - t172 * g(1); t6 * mrSges(5,1) + t2 * mrSges(6,1) - t5 * mrSges(5,2) - t3 * mrSges(6,2) - t15 * t91 - t28 * t92 + t29 * t90 + t216 * t78 + t214 * t77 + (Ifges(6,3) + Ifges(5,3)) * qJDD(4) + (t100 + t153) * g(1) + (t53 + (-g(1) * t139 + t2) * m(6)) * pkin(4) + (-m(6) * (-t13 + t15) + t89) * t16 + ((-t61 / 0.2e1 - t62 / 0.2e1 - t103 / 0.2e1 - t104 / 0.2e1 - t32 * mrSges(5,2) - t27 * mrSges(6,2) + (-Ifges(6,5) / 0.2e1 - Ifges(5,5) / 0.2e1) * qJD(4) + t157) * t139 + (t59 / 0.2e1 + t60 / 0.2e1 - t32 * mrSges(5,1) - t27 * mrSges(6,1) + t221 * t182 + (Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1) * qJD(4) + (-m(6) * t27 - t75) * pkin(4) + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t181 + t156) * t136) * t131 + (-g(2) * t112 + g(3) * t113) * (t136 * (m(6) * pkin(4) + t197) + t139 * t220); (t136 * t89 - t139 * t91) * t131 + (-t131 * t150 - t155 + t8) * m(6) + t30;];
tau = t1;
