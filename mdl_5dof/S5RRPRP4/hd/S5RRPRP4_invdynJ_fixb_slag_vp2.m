% Calculate vector of inverse dynamics joint torques for
% S5RRPRP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4]';
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
% Datum: 2019-12-31 19:53
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPRP4_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRP4_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:52:37
% EndTime: 2019-12-31 19:52:42
% DurationCPUTime: 2.87s
% Computational Cost: add. (1687->274), mult. (2218->333), div. (0->0), fcn. (898->8), ass. (0->136)
t255 = Ifges(5,1) + Ifges(6,1);
t247 = Ifges(6,4) + Ifges(5,5);
t246 = Ifges(6,6) - Ifges(5,6);
t111 = cos(qJ(4));
t108 = sin(qJ(4));
t207 = Ifges(6,5) * t108;
t209 = Ifges(5,4) * t108;
t254 = t255 * t111 + t207 - t209;
t253 = mrSges(3,1) - mrSges(4,2);
t252 = -mrSges(6,1) - mrSges(5,1);
t251 = mrSges(3,2) - mrSges(4,3);
t250 = mrSges(6,3) - mrSges(5,2);
t104 = qJD(1) + qJD(2);
t189 = t104 * t111;
t173 = mrSges(6,2) * t189;
t235 = -mrSges(5,3) * t189 - t252 * qJD(4) - t173;
t190 = t104 * t108;
t174 = mrSges(6,2) * t190;
t68 = qJD(4) * mrSges(6,3) - t174;
t236 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t190 + t68;
t114 = -pkin(2) - pkin(7);
t112 = cos(qJ(2));
t205 = pkin(1) * qJD(1);
t171 = t112 * t205;
t145 = qJD(3) - t171;
t48 = t104 * t114 + t145;
t199 = t111 * t48;
t22 = -qJD(4) * pkin(4) + qJD(5) - t199;
t202 = t108 * t48;
t28 = qJD(4) * qJ(5) + t202;
t131 = t108 * t22 + t111 * t28;
t183 = qJD(4) * t108;
t103 = qJDD(1) + qJDD(2);
t109 = sin(qJ(2));
t204 = pkin(1) * qJD(2);
t170 = t109 * t204;
t218 = pkin(1) * t112;
t59 = -qJD(1) * t170 + qJDD(1) * t218;
t127 = qJDD(3) - t59;
t23 = t103 * t114 + t127;
t6 = t111 * t23 - t183 * t48;
t182 = qJD(4) * t111;
t7 = t108 * t23 + t48 * t182;
t144 = t7 * t108 + t6 * t111;
t4 = qJDD(4) * qJ(5) + qJD(4) * qJD(5) + t7;
t5 = -qJDD(4) * pkin(4) + qJDD(5) - t6;
t241 = t4 * t108 - t111 * t5;
t53 = t103 * t111 - t104 * t183;
t35 = -qJDD(4) * mrSges(6,1) + t53 * mrSges(6,2);
t54 = t103 * t108 + t104 * t182;
t37 = -mrSges(6,2) * t54 + qJDD(4) * mrSges(6,3);
t239 = m(5) * t144 + m(6) * (qJD(4) * t131 + t241) + (-qJDD(4) * mrSges(5,2) - mrSges(5,3) * t54 + t37) * t108 + (qJDD(4) * mrSges(5,1) - mrSges(5,3) * t53 - t35) * t111;
t249 = (-t235 * t108 + t236 * t111) * qJD(4) + t239;
t248 = -mrSges(6,2) - mrSges(5,3) - t253;
t245 = t247 * qJD(4) + t254 * t104;
t208 = Ifges(5,4) * t111;
t244 = t108 * (Ifges(6,3) * t111 - t207) + t111 * (-Ifges(5,1) * t108 - t208);
t139 = t111 * mrSges(6,1) + t108 * mrSges(6,3);
t141 = mrSges(5,1) * t111 - mrSges(5,2) * t108;
t172 = t109 * t205;
t193 = qJ(5) * t111;
t155 = -t108 * pkin(4) + t193;
t69 = qJ(3) - t155;
t21 = t104 * t69 + t172;
t65 = qJ(3) * t104 + t172;
t243 = t21 * t139 + t65 * t141;
t242 = -t247 * t108 + t246 * t111;
t238 = t253 * t104 - t236 * t108 - t235 * t111;
t140 = mrSges(5,1) * t108 + mrSges(5,2) * t111;
t52 = t140 * t104;
t234 = mrSges(4,3) * t104 + t52;
t107 = qJ(1) + qJ(2);
t100 = cos(t107);
t99 = sin(t107);
t233 = -g(1) * t99 + g(2) * t100;
t191 = t100 * t111;
t192 = t100 * t108;
t228 = t252 * t192 + t250 * t191 + t251 * t100 + ((-m(5) - m(6)) * t114 - t248) * t99;
t200 = t108 * t99;
t227 = t248 * t100 + t252 * t200 + (t111 * t250 + t251) * t99;
t224 = -m(4) - m(5);
t223 = pkin(2) * t99;
t221 = t111 / 0.2e1;
t220 = pkin(1) * t109;
t110 = sin(qJ(1));
t219 = pkin(1) * t110;
t113 = cos(qJ(1));
t102 = t113 * pkin(1);
t211 = t100 * pkin(2) + t99 * qJ(3);
t206 = Ifges(6,5) * t111;
t203 = t103 * mrSges(4,2);
t197 = t112 * t65;
t194 = qJ(3) * t103;
t186 = qJD(2) * t112;
t185 = qJD(3) * t104;
t184 = qJD(4) * t104;
t92 = t100 * pkin(7);
t176 = t92 + t211;
t175 = t102 + t211;
t169 = pkin(1) * t186;
t95 = -pkin(2) - t218;
t84 = t100 * qJ(3);
t161 = t84 - t219;
t160 = -t184 / 0.2e1;
t153 = t104 * t171;
t147 = -qJD(5) * t111 + qJD(3);
t146 = (t108 ^ 2 + t111 ^ 2) * t48 * t109;
t60 = (qJD(1) * t186 + qJDD(1) * t109) * pkin(1);
t117 = t60 + t194;
t29 = t117 + t185;
t80 = qJD(3) + t169;
t91 = qJ(3) + t220;
t143 = t29 * t91 + t65 * t80;
t138 = t108 * mrSges(6,1) - t111 * mrSges(6,3);
t135 = -Ifges(5,2) * t108 + t208;
t132 = pkin(4) * t111 + qJ(5) * t108;
t130 = t29 * qJ(3) + t65 * qJD(3);
t128 = pkin(4) * t192 - qJ(5) * t191 + t84;
t47 = pkin(4) * t182 + qJ(5) * t183 + t147;
t124 = t108 * (-Ifges(5,2) * t111 - t209);
t121 = (t108 * t28 - t111 * t22) * t109;
t120 = pkin(4) * t200 - t193 * t99 + t176;
t2 = pkin(4) * t54 - qJ(5) * t53 + t104 * t147 + t117;
t38 = -pkin(2) * t103 + t127;
t79 = Ifges(6,5) * t189;
t43 = Ifges(6,6) * qJD(4) + Ifges(6,3) * t190 + t79;
t44 = Ifges(5,6) * qJD(4) + t104 * t135;
t115 = t108 * (Ifges(6,5) * t53 + Ifges(6,6) * qJDD(4)) / 0.2e1 - t144 * mrSges(5,3) + t2 * t138 - t60 * mrSges(3,2) + t59 * mrSges(3,1) + t38 * mrSges(4,2) - t44 * t182 / 0.2e1 - t108 * (Ifges(5,4) * t53 + Ifges(5,6) * qJDD(4)) / 0.2e1 + t124 * t160 + (t140 + mrSges(4,3)) * t29 + t254 * t53 / 0.2e1 + (t247 * qJDD(4) + t255 * t53) * t221 + (t246 * t108 + t247 * t111) * qJDD(4) / 0.2e1 + t244 * t184 / 0.2e1 - t245 * t183 / 0.2e1 + (t104 * (-Ifges(6,1) * t108 + t206) + t43) * t182 / 0.2e1 + (Ifges(4,1) + Ifges(3,3)) * t103 + (-t182 * t28 - t183 * t22 - t241) * mrSges(6,2) + (t206 / 0.2e1 - t135 / 0.2e1 + (Ifges(6,5) - Ifges(5,4)) * t221 + (Ifges(6,3) + Ifges(5,2) / 0.2e1) * t108) * t54 + (t243 + t242 * qJD(4) / 0.2e1) * qJD(4);
t63 = -pkin(2) * t104 + t145;
t58 = t69 + t220;
t51 = t138 * t104;
t50 = t132 * t104;
t31 = t47 + t169;
t14 = mrSges(5,1) * t54 + mrSges(5,2) * t53;
t13 = mrSges(6,1) * t54 - mrSges(6,3) * t53;
t1 = [m(3) * (t109 * t60 + t112 * t59) * pkin(1) + m(5) * (t146 * t204 + t143) + m(6) * (t121 * t204 + t2 * t58 + t21 * t31) + t115 + t95 * t203 + t91 * t14 + t31 * t51 + t58 * t13 - t104 * mrSges(3,2) * t169 + m(4) * (t38 * t95 + t143) + Ifges(2,3) * qJDD(1) + t234 * t80 + (mrSges(3,1) * t218 - mrSges(3,2) * t220 + t91 * mrSges(4,3)) * t103 + (-m(3) * t102 - mrSges(2,1) * t113 + mrSges(2,2) * t110 - m(6) * (t102 + t120) - m(5) * (t92 + t175) - m(4) * t175 + t227) * g(2) + (-m(6) * (t128 - t219) - m(4) * (t161 - t223) + m(3) * t219 + mrSges(2,1) * t110 + mrSges(2,2) * t113 - m(5) * t161 + t228) * g(1) + (m(4) * t63 - t238) * t170 + (t236 * t182 - t235 * t183 + t239) * (-pkin(7) + t95); -pkin(2) * t203 + t115 + m(5) * t130 + m(4) * (-pkin(2) * t38 + t130) + t69 * t13 + t47 * t51 + qJD(3) * t52 + qJ(3) * t14 + m(6) * (t2 * t69 + t21 * t47) + mrSges(3,2) * t153 + (-t51 - t52) * t171 + (-m(6) * (t112 * t21 + t121) - m(4) * (t109 * t63 + t197) - m(5) * (t146 + t197)) * t205 + (-t153 + t194 + t185) * mrSges(4,3) + (-m(4) * t211 - m(5) * t176 - m(6) * t120 + t227) * g(2) + (-m(4) * (t84 - t223) - m(6) * t128 - m(5) * t84 + t228) * g(1) + t238 * t172 + t249 * t114; t203 + m(4) * t38 + (-m(6) * t21 + t224 * t65 - t234 - t51) * t104 + t233 * (m(6) - t224) + t249; -(-Ifges(6,1) * t190 + t43 + t79) * t189 / 0.2e1 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + (t140 + t138) * g(3) + (-t5 * pkin(4) - g(3) * t155 + t4 * qJ(5) + t28 * qJD(5) - t131 * t48 - t21 * t50) * m(6) + (t132 * m(6) + t139 + t141) * t233 + t28 * t173 + t22 * t174 + qJD(5) * t68 - t50 * t51 - pkin(4) * t35 + qJ(5) * t37 + t4 * mrSges(6,3) - t5 * mrSges(6,1) + t6 * mrSges(5,1) - t7 * mrSges(5,2) + t44 * t189 / 0.2e1 + t235 * t202 - t236 * t199 + t242 * t160 + t245 * t190 / 0.2e1 + t246 * t54 + t247 * t53 + (-t243 + (t124 / 0.2e1 - t244 / 0.2e1) * t104) * t104; t51 * t189 - qJD(4) * t68 + (-g(3) * t108 - t28 * qJD(4) - t111 * t233 + t21 * t189 + t5) * m(6) + t35;];
tau = t1;
