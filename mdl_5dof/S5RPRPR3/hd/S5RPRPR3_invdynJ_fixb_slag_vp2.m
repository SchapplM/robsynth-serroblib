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
% m [6x1]
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:20:36
% EndTime: 2022-01-23 09:20:41
% DurationCPUTime: 2.47s
% Computational Cost: add. (2271->265), mult. (3897->370), div. (0->0), fcn. (2232->14), ass. (0->141)
t127 = qJD(1) + qJD(3);
t131 = cos(pkin(9));
t99 = -t127 * t131 + qJD(5);
t246 = Ifges(6,5) * t99;
t129 = sin(pkin(9));
t137 = cos(qJ(3));
t130 = sin(pkin(8));
t221 = pkin(1) * t130;
t177 = qJD(1) * t221;
t134 = sin(qJ(3));
t132 = cos(pkin(8));
t115 = pkin(1) * t132 + pkin(2);
t97 = t115 * qJD(1);
t200 = t134 * t97;
t66 = t137 * t177 + t200;
t48 = qJ(4) * t127 + t66;
t35 = -t131 * qJD(2) + t129 * t48;
t207 = t129 * t35;
t36 = qJD(2) * t129 + t131 * t48;
t229 = t131 * t36;
t156 = t207 + t229;
t125 = t129 ^ 2;
t243 = mrSges(5,3) * (t131 ^ 2 + t125);
t240 = t127 * t243;
t245 = -m(5) * t156 - t240;
t133 = sin(qJ(5));
t244 = -t133 / 0.2e1;
t242 = -mrSges(5,3) + mrSges(4,2);
t158 = -t131 * mrSges(5,1) + t129 * mrSges(5,2);
t241 = mrSges(4,1) - t158;
t136 = cos(qJ(5));
t213 = Ifges(6,4) * t136;
t145 = (-t133 * Ifges(6,2) + t213) * t129;
t214 = Ifges(6,4) * t133;
t146 = (t136 * Ifges(6,1) - t214) * t129;
t65 = -t134 * t177 + t137 * t97;
t160 = qJD(4) - t65;
t91 = -pkin(4) * t131 - pkin(7) * t129 - pkin(3);
t30 = t91 * t127 + t160;
t8 = -t133 * t36 + t136 * t30;
t9 = t133 * t30 + t136 * t36;
t161 = -t133 * t8 + t136 * t9;
t239 = -t161 * mrSges(6,3) + (t127 * t146 + t246) * t244 + (-t35 * mrSges(6,2) - t246 / 0.2e1) * t133 + (t35 * mrSges(6,1) - Ifges(6,6) * t99 - t127 * t145 / 0.2e1) * t136;
t176 = qJD(3) * t221;
t237 = -qJD(1) * t176 + t115 * qJDD(1);
t124 = qJDD(1) + qJDD(3);
t185 = qJD(4) * t127;
t198 = pkin(1) * qJDD(1);
t169 = t130 * t198;
t186 = qJD(3) * t137;
t28 = t237 * t134 + t137 * t169 + t97 * t186;
t19 = qJ(4) * t124 + t185 + t28;
t13 = -t131 * qJDD(2) + t129 * t19;
t14 = qJDD(2) * t129 + t131 * t19;
t202 = t131 * t14;
t236 = t129 * t13 + t202;
t234 = t136 * (-Ifges(6,1) * t133 - t213) / 0.2e1 + (-Ifges(6,2) * t136 - t214) * t244;
t184 = qJD(4) * t131;
t189 = t131 * t136;
t190 = t131 * t133;
t67 = -qJ(4) * t190 + t136 * t91;
t232 = t67 * qJD(5) - t133 * t66 + t136 * t184 - t65 * t189;
t68 = qJ(4) * t189 + t133 * t91;
t231 = -t68 * qJD(5) - t133 * t184 - t136 * t66 + t65 * t190;
t128 = qJ(1) + pkin(8);
t121 = qJ(3) + t128;
t114 = cos(t121);
t196 = t114 * t131;
t197 = t114 * t129;
t230 = pkin(4) * t196 + pkin(7) * t197;
t113 = sin(t121);
t60 = t113 * t190 + t114 * t136;
t61 = -t113 * t189 + t114 * t133;
t226 = -t61 * mrSges(6,1) - t60 * mrSges(6,2) + t242 * t114 + (-m(6) * t91 + t129 * mrSges(6,3) + t241) * t113;
t62 = t113 * t136 - t114 * t190;
t63 = t113 * t133 + t114 * t189;
t225 = -t114 * mrSges(4,1) - mrSges(5,1) * t196 - t63 * mrSges(6,1) - t62 * mrSges(6,2) + (mrSges(5,2) - mrSges(6,3)) * t197 + t242 * t113;
t183 = qJD(5) * t127;
t69 = (t124 * t136 - t133 * t183) * t129;
t70 = (-t124 * t133 - t136 * t183) * t129;
t27 = -mrSges(6,1) * t70 + mrSges(6,2) * t69;
t224 = t124 * t243 + t129 * t27;
t223 = -m(5) * t160 + (m(5) * pkin(3) + t241) * t127;
t222 = m(3) + m(4);
t135 = sin(qJ(1));
t220 = pkin(1) * t135;
t219 = pkin(3) * t113;
t138 = cos(qJ(1));
t122 = t138 * pkin(1);
t211 = t127 * mrSges(4,2);
t157 = mrSges(6,1) * t133 + mrSges(6,2) * t136;
t193 = t127 * t129;
t64 = t157 * t193;
t206 = t129 * t64;
t205 = t129 * t65;
t204 = t13 * t131;
t80 = t134 * t115 + t137 * t221;
t195 = t124 * t129;
t194 = t124 * t131;
t192 = t129 * t133;
t191 = t129 * t136;
t188 = t114 * pkin(3) + t113 * qJ(4);
t120 = cos(t128);
t187 = pkin(2) * t120 + t122;
t96 = qJDD(5) - t194;
t180 = Ifges(6,5) * t69 + Ifges(6,6) * t70 + Ifges(6,3) * t96;
t179 = mrSges(6,3) * t192;
t178 = mrSges(6,3) * t191;
t79 = t115 * t137 - t134 * t221;
t164 = t187 + t188;
t78 = -mrSges(5,1) * t194 + mrSges(5,2) * t195;
t119 = sin(t128);
t163 = -pkin(2) * t119 - t220;
t162 = -g(1) * t113 + g(2) * t114;
t58 = -mrSges(6,2) * t99 - t127 * t179;
t59 = mrSges(6,1) * t99 - t127 * t178;
t155 = t133 * t59 - t136 * t58;
t29 = -qJD(3) * t200 - t134 * t169 + t237 * t137;
t72 = t115 * t186 - t134 * t176;
t105 = t114 * qJ(4);
t152 = t105 + t163;
t49 = -t79 + t91;
t76 = qJ(4) + t80;
t21 = t133 * t49 + t76 * t189;
t20 = t136 * t49 - t76 * t190;
t71 = qJD(4) + t72;
t151 = (t13 * t76 + t35 * t71) * t129;
t144 = qJDD(4) - t29;
t25 = -pkin(3) * t124 + t144;
t12 = t91 * t124 + t144;
t3 = t8 * qJD(5) + t12 * t133 + t136 * t14;
t4 = -t9 * qJD(5) + t12 * t136 - t133 * t14;
t82 = t157 * t129;
t140 = t96 * (-Ifges(6,3) * t131 + (Ifges(6,5) * t136 - Ifges(6,6) * t133) * t129) / 0.2e1 + Ifges(4,3) * t124 + t13 * t82 + t25 * t158 - t28 * mrSges(4,2) + t29 * mrSges(4,1) + (Ifges(5,4) * t129 + Ifges(5,2) * t131) * t194 + (Ifges(5,1) * t129 + Ifges(5,4) * t131) * t195 + (Ifges(6,1) * t69 + Ifges(6,4) * t70 + Ifges(6,5) * t96) * t191 / 0.2e1 - (Ifges(6,4) * t69 + Ifges(6,2) * t70 + Ifges(6,6) * t96) * t192 / 0.2e1 - t131 * t180 / 0.2e1 + t4 * (-mrSges(6,1) * t131 - t178) + t3 * (mrSges(6,2) * t131 - t179) + t70 * (-Ifges(6,6) * t131 + t145) / 0.2e1 + t69 * (-Ifges(6,5) * t131 + t146) / 0.2e1 + t234 * t125 * t183 + t236 * mrSges(5,3) + t239 * qJD(5) * t129;
t77 = -pkin(3) - t79;
t73 = t80 * qJD(3);
t38 = -mrSges(6,2) * t96 + mrSges(6,3) * t70;
t37 = mrSges(6,1) * t96 - mrSges(6,3) * t69;
t6 = -t21 * qJD(5) + t136 * t73 - t71 * t190;
t5 = t20 * qJD(5) + t133 * t73 + t71 * t189;
t1 = [m(4) * (t28 * t80 + t29 * t79 + t66 * t72) - t72 * t211 + t140 + t77 * t78 + t5 * t58 + t6 * t59 + t20 * t37 + t21 * t38 - 0.2e1 * mrSges(3,2) * t169 + 0.2e1 * t132 * mrSges(3,1) * t198 + m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t151) + m(5) * (t25 * t77 + t151) + (-m(4) * t65 - t223) * t73 + (m(5) * t202 + t224) * t76 + (m(5) * t229 + t206 + t240) * t71 + (t79 * mrSges(4,1) - t80 * mrSges(4,2)) * t124 + (m(3) * (t130 ^ 2 + t132 ^ 2) * pkin(1) ^ 2 + Ifges(3,3) + Ifges(2,3)) * qJDD(1) + (-mrSges(2,1) * t138 + mrSges(2,2) * t135 - m(5) * t164 - m(3) * t122 - mrSges(3,1) * t120 + mrSges(3,2) * t119 - m(4) * t187 - m(6) * (t164 + t230) + t225) * g(2) + (mrSges(2,1) * t135 + mrSges(2,2) * t138 - m(4) * t163 - m(5) * (t152 - t219) + m(3) * t220 + mrSges(3,1) * t119 + mrSges(3,2) * t120 - m(6) * t152 + t226) * g(1); -t131 * t27 + t222 * qJDD(2) + (-t133 * t37 + t136 * t38 + (-t133 * t58 - t136 * t59) * qJD(5)) * t129 + m(5) * (t129 * t14 - t204) + m(6) * (-t204 + (-t133 * t4 + t136 * t3 + (-t133 * t9 - t136 * t8) * qJD(5)) * t129) + (-m(5) - m(6) - t222) * g(3); qJD(4) * t206 - t64 * t205 + m(5) * (-pkin(3) * t25 + t156 * qJD(4)) + t140 + t67 * t37 + t68 * t38 - pkin(3) * t78 + t231 * t59 + t232 * t58 + t185 * t243 + (m(5) * t236 + t224) * qJ(4) + (t3 * t68 + t4 * t67 + (qJ(4) * t13 + qJD(4) * t35) * t129 - t35 * t205 + t232 * t9 + t231 * t8) * m(6) + t223 * t66 + (t211 + t245) * t65 + (-m(6) * (t188 + t230) - m(5) * t188 + t225) * g(2) + (-m(6) * t105 - m(5) * (t105 - t219) + t226) * g(1); t133 * t38 + t136 * t37 - t155 * qJD(5) + (t161 * qJD(5) + t133 * t3 + t136 * t4 + t162) * m(6) + (t162 + t25) * m(5) + t78 + (-t206 + t155 * t131 - m(6) * (t9 * t189 - t8 * t190 + t207) + t245) * t127; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t58 + t9 * t59 - g(1) * (mrSges(6,1) * t62 - mrSges(6,2) * t63) - g(2) * (-mrSges(6,1) * t60 + mrSges(6,2) * t61) + g(3) * t82 + (-t234 * t193 - t239) * t193 + t180;];
tau = t1;
