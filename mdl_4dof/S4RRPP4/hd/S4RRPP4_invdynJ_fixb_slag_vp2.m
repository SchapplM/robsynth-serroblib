% Calculate vector of inverse dynamics joint torques for
% S4RRPP4
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
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tau [4x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:59
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S4RRPP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(4,1),zeros(3,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [4 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: qJDD has to be [4x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP4_invdynJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:58:53
% EndTime: 2019-12-31 16:59:01
% DurationCPUTime: 4.61s
% Computational Cost: add. (647->270), mult. (1502->325), div. (0->0), fcn. (648->4), ass. (0->132)
t187 = Ifges(4,1) + Ifges(5,1);
t185 = Ifges(5,4) + Ifges(4,5);
t186 = Ifges(4,4) + Ifges(3,5);
t171 = -Ifges(5,5) + t186;
t76 = cos(qJ(2));
t143 = t76 * mrSges(4,3);
t144 = t76 * mrSges(5,2);
t74 = sin(qJ(2));
t63 = t74 * qJ(3);
t107 = pkin(1) + t63;
t68 = t76 * pkin(2);
t89 = -t107 - t68;
t25 = t89 * qJD(1);
t136 = qJD(1) * t74;
t60 = pkin(5) * t136;
t38 = -qJD(2) * pkin(2) + qJD(3) + t60;
t135 = qJD(1) * t76;
t61 = pkin(5) * t135;
t73 = qJD(2) * qJ(3);
t43 = t61 + t73;
t160 = pkin(2) + pkin(3);
t6 = qJD(4) + (t160 * t76 + t107) * qJD(1);
t196 = -pkin(5) * (t38 * t76 - t43 * t74) * m(4) - t25 * (t74 * mrSges(4,1) - t143) - t6 * (-t74 * mrSges(5,1) + t144);
t195 = Ifges(5,2) + Ifges(4,3);
t194 = Ifges(3,6) - Ifges(4,6);
t183 = Ifges(4,6) - Ifges(5,6);
t151 = Ifges(4,5) * t76;
t152 = Ifges(5,4) * t76;
t193 = t187 * t74 - t151 - t152;
t75 = sin(qJ(1));
t77 = cos(qJ(1));
t173 = g(1) * t77 + g(2) * t75;
t188 = -m(5) - m(4);
t192 = -m(3) + t188;
t191 = t185 * t136;
t140 = t68 + t63;
t124 = qJD(1) * qJD(2);
t37 = qJDD(1) * t74 + t76 * t124;
t28 = t37 * pkin(5);
t115 = qJDD(3) + t28;
t9 = -qJDD(2) * pkin(2) + t115;
t190 = -t43 * qJD(2) + t9;
t100 = mrSges(3,1) * t76 - t74 * mrSges(3,2);
t98 = t76 * mrSges(4,1) + t74 * mrSges(4,3);
t189 = -t100 - t98;
t182 = qJD(2) * t183 - t135 * t195 + t191;
t117 = mrSges(4,2) * t136;
t181 = mrSges(3,3) * t136 + t117 + (-mrSges(3,1) - mrSges(4,1)) * qJD(2);
t118 = mrSges(5,3) * t135;
t44 = qJD(2) * mrSges(5,2) - t118;
t116 = mrSges(4,2) * t135;
t46 = qJD(2) * mrSges(4,3) + t116;
t180 = t44 + t46;
t179 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t135 + t46;
t64 = t74 * mrSges(5,2);
t141 = t76 * mrSges(5,1) + t64;
t178 = qJDD(2) * qJ(3) + qJD(2) * qJD(3);
t176 = t186 * t76 - t194 * t74;
t175 = t185 * t74;
t114 = t74 * t124;
t125 = t76 * qJDD(1);
t56 = pkin(5) * t125;
t27 = -pkin(5) * t114 + t56;
t174 = t27 * t76 + t28 * t74;
t127 = qJ(4) * qJD(1);
t31 = t127 * t74 - t60;
t172 = -t31 + qJD(3);
t59 = Ifges(3,4) * t135;
t170 = Ifges(3,1) * t136 + t193 * qJD(1) + t171 * qJD(2) + t59;
t169 = -mrSges(2,1) - t64 + t189;
t138 = qJ(3) * t76;
t168 = t173 * t138;
t167 = m(5) * qJ(4) + mrSges(2,2) - mrSges(4,2) - mrSges(3,3);
t159 = pkin(5) * t74;
t158 = pkin(5) * t76;
t7 = t27 + t178;
t155 = t7 * t76;
t154 = Ifges(3,4) * t74;
t153 = Ifges(3,4) * t76;
t148 = t37 * mrSges(5,3);
t142 = pkin(5) - qJ(4);
t134 = qJD(2) * t74;
t133 = qJD(2) * t76;
t132 = qJD(3) * t74;
t131 = qJD(4) * t74;
t130 = qJD(4) * t76;
t129 = qJDD(1) * pkin(1);
t122 = t76 * pkin(3) + t140;
t121 = pkin(5) * t134;
t119 = mrSges(5,3) * t136;
t49 = t142 * t76;
t36 = t114 - t125;
t108 = -t36 * mrSges(5,1) + t37 * mrSges(5,2);
t106 = -t124 / 0.2e1;
t105 = t124 / 0.2e1;
t103 = -m(5) * t160 - mrSges(5,1);
t35 = -t127 * t76 + t61;
t99 = mrSges(3,1) * t74 + mrSges(3,2) * t76;
t95 = t76 * Ifges(3,2) + t154;
t92 = Ifges(5,5) * t76 + Ifges(5,6) * t74;
t91 = pkin(2) * t74 - t138;
t88 = pkin(1) * t99;
t86 = -t160 * t74 + t138;
t84 = t74 * (Ifges(3,1) * t76 - t154);
t83 = t76 * (Ifges(5,2) * t74 + t152);
t82 = t76 * (Ifges(4,3) * t74 + t151);
t79 = qJ(3) * t37 + qJD(1) * t132 + t129;
t47 = t142 * t74;
t40 = -qJD(2) * mrSges(5,1) - t119;
t39 = -pkin(1) - t140;
t34 = t91 * qJD(1);
t33 = t141 * qJD(1);
t32 = t98 * qJD(1);
t30 = t37 * mrSges(4,2);
t26 = pkin(1) + t122;
t21 = Ifges(3,6) * qJD(2) + qJD(1) * t95;
t18 = t35 + t73;
t17 = qJD(2) * t49 - t131;
t16 = qJD(2) * t91 - t132;
t15 = -t134 * t142 - t130;
t14 = -mrSges(4,2) * t36 + qJDD(2) * mrSges(4,3);
t13 = -qJDD(2) * mrSges(4,1) + t30;
t12 = -qJDD(2) * mrSges(5,1) - t148;
t11 = qJDD(2) * mrSges(5,2) + mrSges(5,3) * t36;
t10 = t86 * qJD(1);
t8 = -qJD(2) * t160 + t172;
t5 = qJD(2) * t86 + t132;
t4 = pkin(2) * t36 - t79;
t3 = qJ(4) * t36 + t56 + (-t121 - t130) * qJD(1) + t178;
t2 = -qJ(4) * t37 - qJD(1) * t131 - qJDD(2) * t160 + t115;
t1 = -t160 * t36 + qJDD(4) + t79;
t19 = [(-qJDD(2) * mrSges(3,1) + t13) * t159 + (t181 * pkin(5) - t8 * mrSges(5,3) + t170 / 0.2e1) * t133 + ((Ifges(3,1) + t187) * t37 + (-Ifges(3,4) + t185) * t36 + t171 * qJDD(2)) * t74 / 0.2e1 - t179 * t121 + t175 * t36 / 0.2e1 - (t183 * qJDD(2) + t185 * t37) * t76 / 0.2e1 + m(4) * (t16 * t25 + t39 * t4 + (t9 * t74 + t155) * pkin(5)) + (-t158 * t36 + t159 * t37 + t174) * mrSges(3,3) + m(3) * (qJDD(1) * pkin(1) ^ 2 + pkin(5) * t174) + t100 * t129 + (-t2 * t74 - t3 * t76 + t173) * mrSges(5,3) + ((t187 * t76 + t175) * t74 + t76 * (-Ifges(3,2) * t74 + t153) + t84) * t105 + m(5) * (t1 * t26 + t15 * t18 + t17 * t8 + t2 * t47 + t3 * t49 + t5 * t6) + (-qJDD(2) * mrSges(3,2) + t14) * t158 + (t133 * t38 + t190 * t74 + t155) * mrSges(4,2) + (t74 * Ifges(3,1) + t153 + t193) * t37 / 0.2e1 + (t186 * t74 + t194 * t76) * qJDD(2) / 0.2e1 - qJDD(2) * (t74 * Ifges(5,5) - Ifges(5,6) * t76) / 0.2e1 + t76 * (Ifges(3,4) * t37 - Ifges(3,2) * t36 + Ifges(3,6) * qJDD(2)) / 0.2e1 + ((t176 / 0.2e1 - t92 / 0.2e1) * qJD(2) - t196) * qJD(2) + (t83 + t82) * t106 - pkin(1) * (mrSges(3,1) * t36 + mrSges(3,2) * t37) + t39 * (mrSges(4,1) * t36 - mrSges(4,3) * t37) + t17 * t40 + t15 * t44 + t47 * t12 + t49 * t11 - t16 * t32 + t5 * t33 + (t167 * t75 + t192 * (t77 * pkin(1) + t75 * pkin(5)) + (t188 * t140 - (m(5) * pkin(3) + mrSges(5,1)) * t76 + t169) * t77) * g(2) + ((m(3) * pkin(1) - m(4) * t89 + m(5) * t107 - t103 * t76 - t169) * t75 + (pkin(5) * t192 + t167) * t77) * g(1) + (t18 * mrSges(5,3) + t182 / 0.2e1 - t21 / 0.2e1) * t134 - t36 * t95 / 0.2e1 - t4 * t98 + t26 * t108 + Ifges(2,3) * qJDD(1) - t88 * t124 + t1 * t141 - t195 * t76 * t36; (t11 + t14) * qJ(3) + t180 * qJD(3) - t181 * t61 + (-Ifges(3,6) + t183) * t36 + t176 * t106 + t179 * t60 - (-Ifges(3,2) * t136 + t170 + t59) * t135 / 0.2e1 + t171 * t37 + t43 * t117 + t8 * t118 - t160 * t12 + (Ifges(5,3) + Ifges(4,2) + Ifges(3,3)) * qJDD(2) + (t99 - t144 - t143 + (m(4) * pkin(2) + mrSges(4,1) - t103) * t74) * t173 + (-m(4) * t140 - m(5) * t122 - t141 + t189) * g(3) + (t3 * qJ(3) - t10 * t6 - t160 * t2 + t172 * t18 - t35 * t8 - t168) * m(5) + (-pkin(2) * t9 + qJ(3) * t7 + qJD(3) * t43 - t25 * t34 - t168) * m(4) + ((t88 + t82 / 0.2e1 + t83 / 0.2e1 - t84 / 0.2e1) * qJD(1) + t196) * qJD(1) - t35 * t40 - t31 * t44 - t10 * t33 + t34 * t32 - t27 * mrSges(3,2) - t28 * mrSges(3,1) - t9 * mrSges(4,1) - pkin(2) * t13 + t7 * mrSges(4,3) - t2 * mrSges(5,1) + t3 * mrSges(5,2) - (t187 * t135 + t182 + t191) * t136 / 0.2e1 + t92 * t105 - t38 * t116 - t18 * t119 + t21 * t136 / 0.2e1; -t148 + t30 + (-mrSges(4,1) - mrSges(5,1)) * qJDD(2) - t180 * qJD(2) - t188 * t76 * g(3) + ((-t32 - t33) * qJD(1) + t188 * t173) * t74 + (-qJD(2) * t18 - t136 * t6 + t2) * m(5) + (t136 * t25 + t190) * m(4); (t40 * t74 + t44 * t76) * qJD(1) + (t1 + g(1) * t75 - g(2) * t77 - (-t18 * t76 - t74 * t8) * qJD(1)) * m(5) + t108;];
tau = t19;
