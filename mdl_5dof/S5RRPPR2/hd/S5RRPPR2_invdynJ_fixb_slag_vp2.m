% Calculate vector of inverse dynamics joint torques for
% S5RRPPR2
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2022-01-20 10:06
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5RRPPR2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-20 10:05:27
% EndTime: 2022-01-20 10:05:35
% DurationCPUTime: 2.81s
% Computational Cost: add. (2409->283), mult. (3834->394), div. (0->0), fcn. (2190->14), ass. (0->154)
t141 = sin(qJ(5));
t259 = -t141 / 0.2e1;
t137 = sin(pkin(9));
t133 = t137 ^ 2;
t139 = cos(pkin(9));
t258 = mrSges(5,3) * (t139 ^ 2 + t133);
t135 = qJD(1) + qJD(2);
t138 = sin(pkin(8));
t140 = cos(pkin(8));
t142 = sin(qJ(2));
t226 = pkin(1) * qJD(1);
t189 = t142 * t226;
t145 = cos(qJ(2));
t188 = t145 * t226;
t95 = pkin(2) * t135 + t188;
t53 = t138 * t95 + t140 * t189;
t47 = qJ(4) * t135 + t53;
t37 = -t139 * qJD(3) + t137 * t47;
t218 = t137 * t37;
t38 = qJD(3) * t137 + t139 * t47;
t246 = t139 * t38;
t164 = t218 + t246;
t243 = t135 * t258;
t257 = -m(5) * t164 - t243;
t256 = -mrSges(5,3) + mrSges(4,2);
t166 = -t139 * mrSges(5,1) + t137 * mrSges(5,2);
t255 = mrSges(4,1) - t166;
t103 = -t135 * t139 + qJD(5);
t144 = cos(qJ(5));
t227 = Ifges(6,4) * t144;
t151 = (-t141 * Ifges(6,2) + t227) * t137;
t228 = Ifges(6,4) * t141;
t152 = (t144 * Ifges(6,1) - t228) * t137;
t160 = -pkin(4) * t139 - pkin(7) * t137 - pkin(3);
t104 = t138 * t189;
t52 = t140 * t95 - t104;
t171 = qJD(4) - t52;
t25 = t135 * t160 + t171;
t8 = -t141 * t38 + t144 * t25;
t9 = t141 * t25 + t144 * t38;
t172 = -t141 * t8 + t144 * t9;
t254 = t37 * (mrSges(6,1) * t144 - mrSges(6,2) * t141) - t172 * mrSges(6,3) + (-t144 * t151 / 0.2e1 + t152 * t259) * t135 + (0.2e1 * Ifges(6,5) * t259 - Ifges(6,6) * t144) * t103;
t132 = qJDD(1) + qJDD(2);
t199 = qJD(4) * t135;
t238 = pkin(1) * t142;
t187 = qJD(2) * t238;
t236 = pkin(1) * t145;
t92 = -qJD(1) * t187 + qJDD(1) * t236;
t75 = pkin(2) * t132 + t92;
t200 = qJD(2) * t145;
t93 = (qJD(1) * t200 + qJDD(1) * t142) * pkin(1);
t36 = t138 * t75 + t140 * t93;
t22 = qJ(4) * t132 + t199 + t36;
t16 = -t139 * qJDD(3) + t137 * t22;
t17 = qJDD(3) * t137 + t139 * t22;
t213 = t139 * t17;
t252 = t137 * t16 + t213;
t250 = t144 * (-Ifges(6,1) * t141 - t227) / 0.2e1 + (-Ifges(6,2) * t144 - t228) * t259;
t198 = qJD(4) * t139;
t203 = t139 * t141;
t234 = pkin(2) * t138;
t115 = qJ(4) + t234;
t202 = t139 * t144;
t233 = pkin(2) * t140;
t89 = t160 - t233;
t44 = t115 * t202 + t141 * t89;
t201 = t140 * t142;
t156 = pkin(1) * (t138 * t145 + t201);
t78 = qJD(1) * t156;
t80 = t140 * t188 - t104;
t248 = -qJD(5) * t44 - t141 * t198 - t144 * t78 + t203 * t80;
t43 = -t115 * t203 + t144 * t89;
t247 = qJD(5) * t43 - t141 * t78 + t144 * t198 - t202 * t80;
t244 = (mrSges(3,1) * t238 + mrSges(3,2) * t236) * t135;
t197 = qJD(5) * t135;
t70 = (t132 * t144 - t141 * t197) * t137;
t71 = (-t132 * t141 - t144 * t197) * t137;
t23 = -mrSges(6,1) * t71 + mrSges(6,2) * t70;
t242 = t132 * t258 + t137 * t23;
t136 = qJ(1) + qJ(2);
t127 = pkin(8) + t136;
t117 = sin(t127);
t118 = cos(t127);
t128 = sin(t136);
t129 = cos(t136);
t64 = t117 * t203 + t118 * t144;
t65 = -t117 * t202 + t118 * t141;
t241 = t128 * mrSges(3,1) - t65 * mrSges(6,1) + t129 * mrSges(3,2) - t64 * mrSges(6,2) + t256 * t118 + (-m(6) * t160 + t137 * mrSges(6,3) + t255) * t117;
t209 = t118 * t139;
t210 = t118 * t137;
t66 = t117 * t144 - t118 * t203;
t67 = t117 * t141 + t118 * t202;
t240 = -t129 * mrSges(3,1) - t118 * mrSges(4,1) - mrSges(5,1) * t209 - t67 * mrSges(6,1) + t128 * mrSges(3,2) - t66 * mrSges(6,2) + (mrSges(5,2) - mrSges(6,3)) * t210 + t256 * t117;
t239 = m(4) * t52 - m(5) * t171 + (m(5) * pkin(3) + t255) * t135;
t143 = sin(qJ(1));
t237 = pkin(1) * t143;
t235 = pkin(2) * t128;
t120 = pkin(2) * t129;
t146 = cos(qJ(1));
t130 = t146 * pkin(1);
t225 = t132 * mrSges(4,1);
t224 = t132 * mrSges(4,2);
t222 = t135 * mrSges(4,2);
t165 = mrSges(6,1) * t141 + mrSges(6,2) * t144;
t206 = t135 * t137;
t69 = t165 * t206;
t217 = t137 * t69;
t216 = t137 * t80;
t214 = t139 * t16;
t208 = t132 * t137;
t207 = t132 * t139;
t205 = t137 * t141;
t204 = t137 * t144;
t121 = pkin(2) + t236;
t84 = pkin(1) * t201 + t138 * t121;
t101 = qJDD(5) - t207;
t194 = Ifges(6,5) * t70 + Ifges(6,6) * t71 + Ifges(6,3) * t101;
t191 = mrSges(6,3) * t205;
t190 = mrSges(6,3) * t204;
t182 = t118 * pkin(3) + t117 * qJ(4) + t120;
t108 = t118 * qJ(4);
t179 = t108 - t235;
t35 = -t138 * t93 + t140 * t75;
t83 = t121 * t140 - t138 * t238;
t82 = -mrSges(5,1) * t207 + mrSges(5,2) * t208;
t174 = -t235 - t237;
t173 = -g(1) * t117 + g(2) * t118;
t170 = qJDD(4) - t35;
t169 = pkin(4) * t209 + pkin(7) * t210 + t182;
t62 = -mrSges(6,2) * t103 - t135 * t191;
t63 = mrSges(6,1) * t103 - t135 * t190;
t163 = t141 * t63 - t144 * t62;
t81 = t140 * pkin(1) * t200 - t138 * t187;
t159 = -pkin(3) * t117 + t179;
t51 = t160 - t83;
t76 = qJ(4) + t84;
t21 = t141 * t51 + t202 * t76;
t20 = t144 * t51 - t203 * t76;
t72 = qJD(4) + t81;
t158 = (t16 * t76 + t37 * t72) * t137;
t24 = -pkin(3) * t132 + t170;
t15 = t132 * t160 + t170;
t3 = qJD(5) * t8 + t141 * t15 + t144 * t17;
t4 = -qJD(5) * t9 - t141 * t17 + t144 * t15;
t88 = t165 * t137;
t147 = t24 * t166 + (Ifges(5,4) * t137 + Ifges(5,2) * t139) * t207 + (Ifges(5,1) * t137 + Ifges(5,4) * t139) * t208 + (Ifges(6,1) * t70 + Ifges(6,4) * t71 + Ifges(6,5) * t101) * t204 / 0.2e1 - (Ifges(6,4) * t70 + Ifges(6,2) * t71 + Ifges(6,6) * t101) * t205 / 0.2e1 - t139 * t194 / 0.2e1 + t4 * (-mrSges(6,1) * t139 - t190) + t3 * (mrSges(6,2) * t139 - t191) + t70 * (-Ifges(6,5) * t139 + t152) / 0.2e1 + t71 * (-Ifges(6,6) * t139 + t151) / 0.2e1 + t101 * (-Ifges(6,3) * t139 + (Ifges(6,5) * t144 - Ifges(6,6) * t141) * t137) / 0.2e1 + t16 * t88 + t92 * mrSges(3,1) - t93 * mrSges(3,2) - t36 * mrSges(4,2) + t35 * mrSges(4,1) + t250 * t133 * t197 + (Ifges(3,3) + Ifges(4,3)) * t132 + t252 * mrSges(5,3) + t254 * qJD(5) * t137;
t119 = -pkin(3) - t233;
t79 = qJD(2) * t156;
t77 = -pkin(3) - t83;
t40 = -mrSges(6,2) * t101 + mrSges(6,3) * t71;
t39 = mrSges(6,1) * t101 - mrSges(6,3) * t70;
t6 = -qJD(5) * t21 + t144 * t79 - t203 * t72;
t5 = qJD(5) * t20 + t141 * t79 + t202 * t72;
t1 = [m(6) * (t20 * t4 + t21 * t3 + t5 * t9 + t6 * t8 + t158) + m(5) * (t24 * t77 + t158) - t81 * t222 - t84 * t224 + m(3) * (t142 * t93 + t145 * t92) * pkin(1) + m(4) * (t35 * t83 + t36 * t84 + t53 * t81) + t77 * t82 + t5 * t62 + t6 * t63 + t20 * t39 + t21 * t40 + t147 + t83 * t225 + Ifges(2,3) * qJDD(1) - t239 * t79 + (m(5) * t213 + t242) * t76 + (m(5) * t246 + t217 + t243) * t72 + (mrSges(3,1) * t236 - mrSges(3,2) * t238) * t132 - t244 * qJD(2) + (-m(3) * t130 - m(6) * (t130 + t169) - m(5) * (t130 + t182) - m(4) * (t120 + t130) - t146 * mrSges(2,1) + t143 * mrSges(2,2) + t240) * g(2) + (-m(6) * (t108 + t174) - m(4) * t174 + t143 * mrSges(2,1) + t146 * mrSges(2,2) + m(3) * t237 - m(5) * (t159 - t237) + t241) * g(1); -t224 * t234 - t69 * t216 + m(5) * (t164 * qJD(4) + t119 * t24) + t119 * t82 + m(4) * (t138 * t36 + t140 * t35) * pkin(2) + t43 * t39 + t44 * t40 + t147 + t225 * t233 + qJD(4) * t217 + t248 * t63 + t247 * t62 + t199 * t258 + t244 * qJD(1) + (m(5) * t252 + t242) * t115 + (t3 * t44 + t4 * t43 + (qJD(4) * t37 + t115 * t16) * t137 - t37 * t216 + t247 * t9 + t248 * t8) * m(6) + t239 * t78 + (-m(4) * t53 + t222 + t257) * t80 + (-m(4) * t120 - m(5) * t182 - m(6) * t169 + t240) * g(2) + (m(4) * t235 - m(5) * t159 - m(6) * t179 + t241) * g(1); m(4) * qJDD(3) - t139 * t23 + (-t141 * t39 + t144 * t40 + (-t141 * t62 - t144 * t63) * qJD(5)) * t137 + m(5) * (t137 * t17 - t214) + m(6) * (-t214 + (-t141 * t4 + t144 * t3 + (-t141 * t9 - t144 * t8) * qJD(5)) * t137) + (-m(4) - m(5) - m(6)) * g(3); t141 * t40 + t144 * t39 - t163 * qJD(5) + (t172 * qJD(5) + t141 * t3 + t144 * t4 + t173) * m(6) + (t173 + t24) * m(5) + t82 + (-t217 + t163 * t139 - m(6) * (t202 * t9 - t203 * t8 + t218) + t257) * t135; -t3 * mrSges(6,2) + t4 * mrSges(6,1) - t8 * t62 + t9 * t63 - g(1) * (mrSges(6,1) * t66 - mrSges(6,2) * t67) - g(2) * (-mrSges(6,1) * t64 + mrSges(6,2) * t65) + g(3) * t88 + (-t250 * t206 - t254) * t206 + t194;];
tau = t1;
