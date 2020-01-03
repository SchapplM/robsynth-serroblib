% Calculate vector of inverse dynamics joint torques for
% S5RRPRP2
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
% Datum: 2019-12-31 19:50
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:49:26
% EndTime: 2019-12-31 19:49:32
% DurationCPUTime: 2.89s
% Computational Cost: add. (2039->291), mult. (3296->362), div. (0->0), fcn. (1657->12), ass. (0->139)
t266 = -mrSges(5,1) - mrSges(6,1);
t262 = Ifges(6,4) + Ifges(5,5);
t261 = Ifges(6,6) - Ifges(5,6);
t134 = sin(qJ(4));
t137 = cos(qJ(4));
t130 = qJD(1) + qJD(2);
t201 = t130 * t137;
t186 = mrSges(5,3) * t201;
t190 = mrSges(6,2) * t201;
t88 = qJD(4) * mrSges(6,3) + t190;
t250 = -qJD(4) * mrSges(5,2) + t186 + t88;
t202 = t130 * t134;
t187 = mrSges(5,3) * t202;
t191 = mrSges(6,2) * t202;
t251 = t266 * qJD(4) + t187 + t191;
t265 = t134 * t251 + t137 * t250;
t164 = t137 * mrSges(6,1) + t134 * mrSges(6,3);
t96 = -t137 * mrSges(5,1) + mrSges(5,2) * t134;
t264 = -t164 + t96;
t263 = -mrSges(6,2) - mrSges(5,3) + mrSges(4,2);
t100 = Ifges(5,4) * t201;
t220 = Ifges(6,5) * t137;
t162 = Ifges(6,1) * t134 - t220;
t260 = Ifges(5,1) * t202 + t262 * qJD(4) + t130 * t162 + t100;
t158 = pkin(4) * t137 + qJ(5) * t134;
t150 = -pkin(3) - t158;
t133 = cos(pkin(8));
t138 = cos(qJ(2));
t198 = qJD(1) * t138;
t188 = pkin(1) * t198;
t89 = pkin(2) * t130 + t188;
t132 = sin(pkin(8));
t135 = sin(qJ(2));
t236 = pkin(1) * t135;
t189 = qJD(1) * t236;
t98 = t132 * t189;
t43 = t133 * t89 - t98;
t16 = t130 * t150 - t43;
t163 = t134 * mrSges(6,1) - t137 * mrSges(6,3);
t165 = mrSges(5,1) * t134 + mrSges(5,2) * t137;
t32 = -pkin(3) * t130 - t43;
t259 = t16 * t163 + t32 * t165;
t258 = t261 * t134 + t262 * t137;
t129 = qJDD(1) + qJDD(2);
t234 = pkin(1) * t138;
t81 = -qJD(2) * t189 + qJDD(1) * t234;
t58 = pkin(2) * t129 + t81;
t82 = (qJD(2) * t198 + qJDD(1) * t135) * pkin(1);
t26 = t132 * t58 + t133 * t82;
t15 = pkin(7) * t129 + t26;
t195 = qJD(4) * t137;
t185 = qJD(3) * t195 + t134 * qJDD(3) + t137 * t15;
t196 = qJD(4) * t134;
t44 = t132 * t89 + t133 * t189;
t33 = pkin(7) * t130 + t44;
t7 = -t196 * t33 + t185;
t28 = qJD(3) * t134 + t137 * t33;
t8 = -qJD(4) * t28 + qJDD(3) * t137 - t134 * t15;
t257 = -t8 * t134 + t137 * t7;
t215 = t134 * t33;
t3 = qJDD(4) * qJ(5) + (qJD(5) - t215) * qJD(4) + t185;
t5 = -qJDD(4) * pkin(4) + qJDD(5) - t8;
t256 = t134 * t5 + t137 * t3;
t25 = -t132 * t82 + t133 * t58;
t14 = -pkin(3) * t129 - t25;
t74 = -t137 * t129 + t130 * t196;
t75 = t129 * t134 + t130 * t195;
t255 = m(5) * t14 + mrSges(5,1) * t74 + mrSges(5,2) * t75;
t27 = qJD(3) * t137 - t215;
t21 = -qJD(4) * pkin(4) + qJD(5) - t27;
t24 = qJD(4) * qJ(5) + t28;
t48 = -qJDD(4) * mrSges(6,1) + t75 * mrSges(6,2);
t252 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t75 - t48;
t49 = -mrSges(6,2) * t74 + qJDD(4) * mrSges(6,3);
t253 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t74 + t49;
t254 = m(5) * ((-t134 * t28 - t137 * t27) * qJD(4) + t257) + m(6) * ((-t134 * t24 + t137 * t21) * qJD(4) + t256) - t196 * t250 + t195 * t251 + t137 * t253 - t134 * t252;
t249 = (mrSges(3,1) * t236 + mrSges(3,2) * t234) * t130;
t131 = qJ(1) + qJ(2);
t125 = pkin(8) + t131;
t113 = sin(t125);
t114 = cos(t125);
t248 = g(1) * t114 + g(2) * t113;
t126 = sin(t131);
t127 = cos(t131);
t243 = mrSges(3,1) * t126 + mrSges(3,2) * t127 + t263 * t114 + (-m(6) * t150 + mrSges(4,1) - t264) * t113;
t205 = t114 * t137;
t206 = t114 * t134;
t242 = -t127 * mrSges(3,1) - t114 * mrSges(4,1) + t126 * mrSges(3,2) + (mrSges(5,2) - mrSges(6,3)) * t206 + t266 * t205 + t263 * t113;
t241 = -m(4) * t43 + m(5) * t32 + (-mrSges(4,1) + t96) * t130;
t155 = -t134 * t27 + t137 * t28;
t156 = t134 * t21 + t137 * t24;
t240 = m(4) * t44 + m(5) * t155 + m(6) * t156 - t130 * mrSges(4,2) + t265;
t237 = t134 / 0.2e1;
t136 = sin(qJ(1));
t235 = pkin(1) * t136;
t233 = pkin(2) * t126;
t118 = pkin(2) * t127;
t232 = pkin(2) * t132;
t231 = pkin(2) * t133;
t139 = cos(qJ(1));
t128 = t139 * pkin(1);
t223 = Ifges(5,4) * t134;
t222 = Ifges(5,4) * t137;
t221 = Ifges(6,5) * t134;
t219 = t129 * mrSges(4,1);
t218 = t129 * mrSges(4,2);
t119 = pkin(2) + t234;
t199 = t133 * t135;
t67 = pkin(1) * t199 + t132 * t119;
t200 = t132 * t135;
t197 = qJD(4) * t130;
t194 = qJD(5) * t134;
t182 = t114 * pkin(3) + t113 * pkin(7) + t118;
t177 = -t197 / 0.2e1;
t107 = t114 * pkin(7);
t174 = t107 - t233;
t66 = -pkin(1) * t200 + t119 * t133;
t169 = -t233 - t235;
t168 = pkin(4) * t205 + qJ(5) * t206 + t182;
t161 = Ifges(5,2) * t137 + t223;
t157 = pkin(4) * t134 - qJ(5) * t137;
t151 = -pkin(3) * t113 + t174;
t147 = pkin(1) * (t132 * t138 + t199);
t146 = t134 * (Ifges(5,1) * t137 - t223);
t145 = t137 * (Ifges(6,3) * t134 + t220);
t70 = pkin(4) * t196 - qJ(5) * t195 - t194;
t62 = qJD(2) * t147;
t99 = Ifges(6,5) * t202;
t54 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t201 + t99;
t55 = Ifges(5,6) * qJD(4) + t130 * t161;
t9 = pkin(4) * t74 - qJ(5) * t75 - t130 * t194 + t14;
t140 = t137 * (Ifges(5,4) * t75 + Ifges(5,6) * qJDD(4)) / 0.2e1 - t9 * t164 - t137 * (Ifges(6,5) * t75 + Ifges(6,6) * qJDD(4)) / 0.2e1 - t55 * t196 / 0.2e1 + t146 * t197 / 0.2e1 + t81 * mrSges(3,1) - t82 * mrSges(3,2) + t14 * t96 + t25 * mrSges(4,1) - t26 * mrSges(4,2) + t145 * t177 + (Ifges(5,1) * t134 + t162 + t222) * t75 / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t75 + t262 * qJDD(4)) * t237 + (t262 * t134 - t261 * t137) * qJDD(4) / 0.2e1 + (t130 * (Ifges(6,1) * t137 + t221) + t54) * t196 / 0.2e1 + (Ifges(3,3) + Ifges(4,3)) * t129 + (-t195 * t27 - t196 * t28 + t257) * mrSges(5,3) + (-t161 / 0.2e1 + t221 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t237 + (-Ifges(5,2) / 0.2e1 - Ifges(6,3)) * t137) * t74 + (t130 * (-Ifges(5,2) * t134 + t222) + t260) * t195 / 0.2e1 + (t21 * t195 - t24 * t196 + t256) * mrSges(6,2) + (t259 + t258 * qJD(4) / 0.2e1) * qJD(4);
t80 = t150 - t231;
t73 = t157 * t130;
t71 = t164 * t130;
t38 = t150 - t66;
t31 = t62 + t70;
t29 = mrSges(6,1) * t74 - mrSges(6,3) * t75;
t1 = [m(4) * (t25 * t66 + t26 * t67) + t66 * t219 - t67 * t218 + m(6) * (t16 * t31 + t38 * t9) + t140 - t31 * t71 + t38 * t29 + Ifges(2,3) * qJDD(1) + t241 * t62 + (mrSges(3,1) * t234 - mrSges(3,2) * t236) * t129 - t249 * qJD(2) + (-m(6) * (t128 + t168) - m(5) * (t128 + t182) - m(4) * (t118 + t128) - mrSges(2,1) * t139 + t136 * mrSges(2,2) - m(3) * t128 + t242) * g(2) + (-m(6) * (t107 + t169) - m(4) * t169 + t136 * mrSges(2,1) + mrSges(2,2) * t139 + m(3) * t235 - m(5) * (t151 - t235) + t243) * g(1) + t255 * (-pkin(3) - t66) + (m(3) * (t135 * t82 + t138 * t81) + t240 * (t133 * t138 - t200) * qJD(2)) * pkin(1) + t254 * (pkin(7) + t67); t219 * t231 + m(4) * (t132 * t26 + t133 * t25) * pkin(2) - t218 * t232 + t140 - t70 * t71 + t80 * t29 + m(6) * (t16 * t70 + t80 * t9) + (-m(4) * t118 - m(5) * t182 - m(6) * t168 + t242) * g(2) + (m(4) * t233 - m(5) * t151 - m(6) * t174 + t243) * g(1) - t240 * (t133 * t188 - t98) + t255 * (-pkin(3) - t231) + (t249 + (-m(6) * t16 - t241 + t71) * t147) * qJD(1) + t254 * (pkin(7) + t232); m(4) * qJDD(3) + t252 * t137 + t253 * t134 + t265 * qJD(4) + m(5) * (qJD(4) * t155 + t134 * t7 + t137 * t8) + m(6) * (qJD(4) * t156 + t134 * t3 - t137 * t5) + (-m(4) - m(5) - m(6)) * g(3); -(-Ifges(5,2) * t202 + t100 + t260) * t201 / 0.2e1 - (Ifges(6,1) * t201 + t54 + t99) * t202 / 0.2e1 + t258 * t177 + t24 * t191 + (-t5 * pkin(4) - g(3) * t158 + t3 * qJ(5) + t24 * qJD(5) - t16 * t73) * m(6) + t262 * t75 + t261 * t74 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) - t21 * t190 + (-m(6) * t24 + t186 - t250) * t27 + (-m(6) * t21 + t187 - t251) * t28 + t55 * t202 / 0.2e1 + qJD(5) * t88 + t73 * t71 - pkin(4) * t48 + qJ(5) * t49 + t3 * mrSges(6,3) - t5 * mrSges(6,1) - t7 * mrSges(5,2) + t8 * mrSges(5,1) + t264 * g(3) + (m(6) * t157 + t163 + t165) * t248 + ((-t146 / 0.2e1 + t145 / 0.2e1) * t130 - t259) * t130; -t71 * t202 - qJD(4) * t88 + (g(3) * t137 - t24 * qJD(4) - t134 * t248 + t16 * t202 + t5) * m(6) + t48;];
tau = t1;
