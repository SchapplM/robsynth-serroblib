% Calculate vector of inverse dynamics joint torques for
% S5RPRPR3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:52
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RPRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:51:14
% EndTime: 2019-12-05 17:51:21
% DurationCPUTime: 2.56s
% Computational Cost: add. (2271->260), mult. (3897->364), div. (0->0), fcn. (2232->14), ass. (0->140)
t119 = qJD(1) + qJD(3);
t123 = cos(pkin(9));
t96 = -t119 * t123 + qJD(5);
t240 = Ifges(6,5) * t96;
t121 = sin(pkin(9));
t129 = cos(qJ(3));
t122 = sin(pkin(8));
t214 = pkin(1) * t122;
t172 = qJD(1) * t214;
t126 = sin(qJ(3));
t124 = cos(pkin(8));
t108 = pkin(1) * t124 + pkin(2);
t94 = t108 * qJD(1);
t191 = t126 * t94;
t66 = t129 * t172 + t191;
t48 = qJ(4) * t119 + t66;
t35 = -qJD(2) * t123 + t121 * t48;
t198 = t121 * t35;
t36 = qJD(2) * t121 + t123 * t48;
t222 = t123 * t36;
t151 = t198 + t222;
t117 = t121 ^ 2;
t237 = mrSges(5,3) * (t123 ^ 2 + t117);
t234 = t119 * t237;
t239 = -m(5) * t151 - t234;
t125 = sin(qJ(5));
t238 = -t125 / 0.2e1;
t228 = t123 * mrSges(5,1) - mrSges(5,2) * t121;
t236 = mrSges(4,1) + t228;
t235 = mrSges(4,2) - mrSges(5,3);
t90 = -pkin(4) * t123 - pkin(7) * t121 - pkin(3);
t233 = -m(6) * t90 + t121 * mrSges(6,3) + t236;
t128 = cos(qJ(5));
t203 = Ifges(6,4) * t128;
t140 = (-t125 * Ifges(6,2) + t203) * t121;
t204 = Ifges(6,4) * t125;
t141 = (t128 * Ifges(6,1) - t204) * t121;
t65 = -t126 * t172 + t129 * t94;
t159 = qJD(4) - t65;
t30 = t119 * t90 + t159;
t8 = -t125 * t36 + t128 * t30;
t9 = t125 * t30 + t128 * t36;
t155 = -t125 * t8 + t128 * t9;
t232 = -t155 * mrSges(6,3) + (t119 * t141 + t240) * t238 + (-t240 / 0.2e1 - t35 * mrSges(6,2)) * t125 + (-Ifges(6,6) * t96 + t35 * mrSges(6,1) - t119 * t140 / 0.2e1) * t128;
t171 = qJD(3) * t214;
t230 = -qJD(1) * t171 + qJDD(1) * t108;
t116 = qJDD(1) + qJDD(3);
t179 = qJD(4) * t119;
t188 = qJDD(1) * pkin(1);
t164 = t122 * t188;
t180 = qJD(3) * t129;
t28 = t126 * t230 + t129 * t164 + t94 * t180;
t19 = qJ(4) * t116 + t179 + t28;
t13 = -qJDD(2) * t123 + t121 * t19;
t14 = qJDD(2) * t121 + t123 * t19;
t193 = t123 * t14;
t229 = t121 * t13 + t193;
t226 = t128 * (-Ifges(6,1) * t125 - t203) / 0.2e1 + (-Ifges(6,2) * t128 - t204) * t238;
t178 = qJD(4) * t123;
t181 = t123 * t128;
t182 = t123 * t125;
t67 = -qJ(4) * t182 + t128 * t90;
t224 = qJD(5) * t67 - t125 * t66 + t128 * t178 - t181 * t65;
t68 = qJ(4) * t181 + t125 * t90;
t223 = -qJD(5) * t68 - t125 * t178 - t128 * t66 + t182 * t65;
t120 = qJ(1) + pkin(8);
t114 = qJ(3) + t120;
t106 = sin(t114);
t107 = cos(t114);
t62 = -t106 * t128 + t107 * t182;
t63 = -t106 * t125 - t107 * t181;
t219 = -t63 * mrSges(6,1) - t62 * mrSges(6,2) - t106 * t235 + t107 * t233;
t60 = t106 * t182 + t107 * t128;
t61 = t106 * t181 - t107 * t125;
t218 = t61 * mrSges(6,1) - t60 * mrSges(6,2) + t106 * t233 + t107 * t235;
t177 = qJD(5) * t119;
t69 = (t116 * t128 - t125 * t177) * t121;
t70 = (-t116 * t125 - t128 * t177) * t121;
t27 = -mrSges(6,1) * t70 + mrSges(6,2) * t69;
t217 = t116 * t237 + t121 * t27;
t216 = -m(5) * t159 + (m(5) * pkin(3) + t236) * t119;
t215 = m(3) + m(4);
t127 = sin(qJ(1));
t213 = pkin(1) * t127;
t130 = cos(qJ(1));
t212 = pkin(1) * t130;
t211 = pkin(3) * t106;
t210 = pkin(3) * t107;
t201 = t119 * mrSges(4,2);
t152 = mrSges(6,1) * t125 + mrSges(6,2) * t128;
t185 = t119 * t121;
t64 = t152 * t185;
t197 = t121 * t64;
t196 = t121 * t65;
t194 = t123 * t13;
t80 = t108 * t126 + t129 * t214;
t189 = qJ(4) * t106;
t187 = t116 * t121;
t186 = t116 * t123;
t184 = t121 * t125;
t183 = t121 * t128;
t93 = qJDD(5) - t186;
t175 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t93;
t174 = mrSges(6,3) * t184;
t173 = mrSges(6,3) * t183;
t79 = t108 * t129 - t126 * t214;
t78 = -mrSges(5,1) * t186 + mrSges(5,2) * t187;
t112 = sin(t120);
t158 = -pkin(2) * t112 - t213;
t113 = cos(t120);
t157 = -pkin(2) * t113 - t212;
t156 = -g(2) * t107 - g(3) * t106;
t58 = -mrSges(6,2) * t96 - t119 * t174;
t59 = mrSges(6,1) * t96 - t119 * t173;
t150 = t125 * t59 - t128 * t58;
t29 = -qJD(3) * t191 - t126 * t164 + t129 * t230;
t72 = t108 * t180 - t126 * t171;
t101 = t107 * qJ(4);
t148 = t101 + t158;
t49 = -t79 + t90;
t76 = qJ(4) + t80;
t21 = t125 * t49 + t181 * t76;
t20 = t128 * t49 - t182 * t76;
t71 = qJD(4) + t72;
t147 = (t13 * t76 + t35 * t71) * t121;
t144 = qJDD(4) - t29;
t138 = t157 - t189;
t25 = -pkin(3) * t116 + t144;
t12 = t116 * t90 + t144;
t3 = qJD(5) * t8 + t12 * t125 + t128 * t14;
t4 = -qJD(5) * t9 + t12 * t128 - t125 * t14;
t82 = t152 * t121;
t132 = (Ifges(5,4) * t121 + Ifges(5,2) * t123) * t186 + (Ifges(5,1) * t121 + Ifges(5,4) * t123) * t187 + (Ifges(6,1) * t69 + Ifges(6,4) * t70 + Ifges(6,5) * t93) * t183 / 0.2e1 - (Ifges(6,4) * t69 + Ifges(6,2) * t70 + Ifges(6,6) * t93) * t184 / 0.2e1 - t123 * t175 / 0.2e1 + t4 * (-mrSges(6,1) * t123 - t173) + t3 * (mrSges(6,2) * t123 - t174) - t25 * t228 + t70 * (-Ifges(6,6) * t123 + t140) / 0.2e1 + t69 * (-Ifges(6,5) * t123 + t141) / 0.2e1 + t93 * (-Ifges(6,3) * t123 + (Ifges(6,5) * t128 - Ifges(6,6) * t125) * t121) / 0.2e1 + Ifges(4,3) * t116 + t13 * t82 + t29 * mrSges(4,1) - t28 * mrSges(4,2) + t226 * t117 * t177 + t229 * mrSges(5,3) + t232 * qJD(5) * t121;
t77 = -pkin(3) - t79;
t73 = t80 * qJD(3);
t38 = -mrSges(6,2) * t93 + mrSges(6,3) * t70;
t37 = mrSges(6,1) * t93 - mrSges(6,3) * t69;
t6 = -qJD(5) * t21 + t128 * t73 - t182 * t71;
t5 = qJD(5) * t20 + t125 * t73 + t181 * t71;
t1 = [-t72 * t201 - 0.2e1 * mrSges(3,2) * t164 + m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t147) + m(5) * (t25 * t77 + t147) + t77 * t78 + t5 * t58 + t6 * t59 + t20 * t37 + t21 * t38 + 0.2e1 * t124 * mrSges(3,1) * t188 + m(4) * (t28 * t80 + t29 * t79 + t66 * t72) + t132 + (-m(4) * t65 - t216) * t73 + (m(5) * t193 + t217) * t76 + (m(5) * t222 + t197 + t234) * t71 + (mrSges(4,1) * t79 - mrSges(4,2) * t80) * t116 + (m(3) * (t122 ^ 2 + t124 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + (-m(4) * t158 - m(6) * t148 + mrSges(2,1) * t127 + mrSges(2,2) * t130 - m(5) * (t148 - t211) + m(3) * t213 + mrSges(3,1) * t112 + mrSges(3,2) * t113 + t218) * g(3) + (-m(5) * (t138 - t210) - m(4) * t157 - m(6) * t138 + mrSges(2,1) * t130 - mrSges(2,2) * t127 + m(3) * t212 + mrSges(3,1) * t113 - mrSges(3,2) * t112 + t219) * g(2); -t123 * t27 + t215 * qJDD(2) + (-t125 * t37 + t128 * t38 + (-t125 * t58 - t128 * t59) * qJD(5)) * t121 + m(5) * (t121 * t14 - t194) + m(6) * (-t194 + (-t125 * t4 + t128 * t3 + (-t125 * t9 - t128 * t8) * qJD(5)) * t121) + (-m(5) - m(6) - t215) * g(1); m(5) * (-pkin(3) * t25 + t151 * qJD(4)) - t64 * t196 - pkin(3) * t78 + t67 * t37 + t68 * t38 + qJD(4) * t197 + t132 + t223 * t59 + t224 * t58 + t179 * t237 + (m(5) * t229 + t217) * qJ(4) + (-t35 * t196 + t3 * t68 + t4 * t67 + (qJ(4) * t13 + qJD(4) * t35) * t121 + t224 * t9 + t223 * t8) * m(6) + t216 * t66 + (t201 + t239) * t65 + (-m(5) * (t101 - t211) - m(6) * t101 + t218) * g(3) + (-m(5) * (-t189 - t210) + m(6) * t189 + t219) * g(2); t125 * t38 + t128 * t37 - t150 * qJD(5) + (t155 * qJD(5) + t125 * t3 + t128 * t4 + t156) * m(6) + (t156 + t25) * m(5) + t78 + (-t197 + t150 * t123 - m(6) * (t181 * t9 - t182 * t8 + t198) + t239) * t119; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t58 + t9 * t59 + g(1) * t82 - g(2) * (mrSges(6,1) * t60 + mrSges(6,2) * t61) - g(3) * (-mrSges(6,1) * t62 + mrSges(6,2) * t63) + (-t226 * t185 - t232) * t185 + t175;];
tau = t1;
