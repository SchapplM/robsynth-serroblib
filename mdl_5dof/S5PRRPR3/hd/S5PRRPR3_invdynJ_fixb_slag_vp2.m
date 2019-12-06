% Calculate vector of inverse dynamics joint torques for
% S5PRRPR3
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
%   pkin=[a2,a3,a4,a5,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:20
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PRRPR3_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR3_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:19:03
% EndTime: 2019-12-05 16:19:15
% DurationCPUTime: 5.13s
% Computational Cost: add. (2946->379), mult. (6552->514), div. (0->0), fcn. (4538->12), ass. (0->172)
t141 = qJDD(3) + qJDD(5);
t144 = qJD(3) + qJD(5);
t149 = sin(qJ(5));
t151 = cos(qJ(5));
t146 = sin(pkin(9));
t147 = cos(pkin(9));
t150 = sin(qJ(3));
t152 = cos(qJ(3));
t104 = -t146 * t150 + t147 * t152;
t98 = t104 * qJD(2);
t184 = t152 * qJD(2);
t185 = t150 * qJD(2);
t99 = -t146 * t184 - t147 * t185;
t169 = t149 * t99 + t151 * t98;
t182 = qJD(2) * qJD(3);
t171 = t150 * t182;
t181 = qJDD(2) * t152;
t110 = -t171 + t181;
t172 = t152 * t182;
t111 = t150 * qJDD(2) + t172;
t63 = t110 * t147 - t111 * t146;
t64 = t110 * t146 + t111 * t147;
t17 = qJD(5) * t169 + t149 * t63 + t151 * t64;
t56 = t149 * t98 - t151 * t99;
t18 = -qJD(5) * t56 - t149 * t64 + t151 * t63;
t138 = t152 * qJDD(1);
t183 = qJD(1) * qJD(3);
t38 = -pkin(6) * t172 + qJDD(3) * pkin(3) - t111 * qJ(4) + t138 + (-pkin(6) * qJDD(2) - qJD(2) * qJD(4) - t183) * t150;
t177 = pkin(6) * t181 + t150 * qJDD(1) + t152 * t183;
t189 = qJD(3) * t150;
t179 = pkin(6) * t189;
t187 = qJD(4) * t152;
t43 = t110 * qJ(4) + (-t179 + t187) * qJD(2) + t177;
t19 = -t146 * t43 + t147 * t38;
t6 = qJDD(3) * pkin(4) - pkin(7) * t64 + t19;
t212 = pkin(7) * t99;
t148 = -qJ(4) - pkin(6);
t120 = t148 * t152;
t186 = t150 * qJD(1);
t95 = -qJD(2) * t120 + t186;
t77 = t146 * t95;
t118 = t148 * t150;
t140 = t152 * qJD(1);
t94 = qJD(2) * t118 + t140;
t86 = qJD(3) * pkin(3) + t94;
t41 = t147 * t86 - t77;
t28 = qJD(3) * pkin(4) + t212 + t41;
t213 = pkin(7) * t98;
t196 = t147 * t95;
t42 = t146 * t86 + t196;
t29 = t42 + t213;
t7 = -t149 * t29 + t151 * t28;
t20 = t146 * t38 + t147 * t43;
t9 = pkin(7) * t63 + t20;
t2 = qJD(5) * t7 + t149 * t6 + t151 * t9;
t207 = Ifges(6,4) * t56;
t50 = Ifges(6,4) * t169;
t24 = Ifges(6,1) * t56 + Ifges(6,5) * t144 + t50;
t8 = t149 * t28 + t151 * t29;
t3 = -qJD(5) * t8 - t149 * t9 + t151 * t6;
t204 = pkin(3) * t152;
t129 = pkin(2) + t204;
t115 = -qJD(2) * t129 + qJD(4);
t65 = -t98 * pkin(4) + t115;
t235 = t3 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,5) * t17 + Ifges(6,6) * t18 + Ifges(6,3) * t141 - (Ifges(6,5) * t169 - Ifges(6,6) * t56) * t144 / 0.2e1 + (t169 * t7 + t56 * t8) * mrSges(6,3) - (-Ifges(6,2) * t56 + t24 + t50) * t169 / 0.2e1 - t65 * (mrSges(6,1) * t56 + mrSges(6,2) * t169) - (Ifges(6,1) * t169 - t207) * t56 / 0.2e1;
t44 = -t146 * t94 - t196;
t30 = t44 - t213;
t46 = t147 * t94 - t77;
t32 = t46 + t212;
t127 = pkin(3) * t147 + pkin(4);
t206 = pkin(3) * t146;
t92 = t127 * t151 - t149 * t206;
t234 = t92 * qJD(5) - t149 * t30 - t151 * t32;
t93 = t127 * t149 + t151 * t206;
t233 = -t93 * qJD(5) + t149 * t32 - t151 * t30;
t143 = pkin(8) + qJ(2);
t133 = sin(t143);
t135 = cos(t143);
t232 = g(1) * t135 + g(2) * t133;
t119 = -t152 * mrSges(4,1) + t150 * mrSges(4,2);
t145 = qJ(3) + pkin(9);
t134 = sin(t145);
t136 = cos(t145);
t139 = qJ(5) + t145;
t125 = sin(t139);
t126 = cos(t139);
t166 = t126 * mrSges(6,1) - t125 * mrSges(6,2);
t231 = t136 * mrSges(5,1) - t134 * mrSges(5,2) - t119 + t166;
t23 = Ifges(6,2) * t169 + Ifges(6,6) * t144 + t207;
t229 = t23 / 0.2e1;
t191 = qJDD(3) / 0.2e1;
t190 = qJDD(2) * pkin(2);
t68 = -pkin(6) * t171 + t177;
t69 = -pkin(6) * t111 - t150 * t183 + t138;
t224 = -t69 * t150 + t68 * t152;
t223 = 0.2e1 * t191;
t170 = pkin(4) * t136 + t204;
t222 = mrSges(3,1) + m(6) * (pkin(2) + t170) + m(5) * t129 + m(4) * pkin(2) + t231;
t221 = mrSges(3,2) + m(6) * (-pkin(7) + t148) - mrSges(6,3) + m(5) * t148 - mrSges(5,3) - m(4) * pkin(6) - mrSges(4,3);
t220 = m(5) * pkin(3);
t217 = t56 / 0.2e1;
t215 = -t99 / 0.2e1;
t214 = m(2) + m(3);
t210 = t150 / 0.2e1;
t208 = Ifges(5,4) * t99;
t205 = pkin(3) * t150;
t201 = -qJD(3) / 0.2e1;
t168 = qJD(3) * t148;
t96 = t150 * t168 + t187;
t97 = -t150 * qJD(4) + t152 * t168;
t47 = t146 * t97 + t147 * t96;
t200 = mrSges(4,2) * t152;
t199 = Ifges(4,4) * t150;
t198 = Ifges(4,4) * t152;
t194 = t152 * Ifges(4,2);
t67 = t146 * t118 - t147 * t120;
t188 = qJD(3) * t152;
t180 = pkin(3) * t189;
t174 = -t63 * mrSges(5,1) + t64 * mrSges(5,2);
t173 = -t18 * mrSges(6,1) + t17 * mrSges(6,2);
t45 = -t146 * t96 + t147 * t97;
t66 = t147 * t118 + t120 * t146;
t165 = -g(1) * t133 + g(2) * t135;
t163 = mrSges(6,1) * t125 + mrSges(6,2) * t126;
t162 = t194 + t199;
t161 = Ifges(4,5) * t152 - Ifges(4,6) * t150;
t105 = t146 * t152 + t147 * t150;
t48 = -pkin(7) * t105 + t66;
t49 = pkin(7) * t104 + t67;
t21 = -t149 * t49 + t151 * t48;
t22 = t149 * t48 + t151 * t49;
t60 = t104 * t151 - t105 * t149;
t61 = t104 * t149 + t105 * t151;
t113 = -pkin(6) * t185 + t140;
t114 = pkin(6) * t184 + t186;
t160 = t113 * t152 + t114 * t150;
t158 = pkin(2) * (mrSges(4,1) * t150 + t200);
t88 = -pkin(3) * t110 + qJDD(4) - t190;
t157 = t150 * (Ifges(4,1) * t152 - t199);
t131 = Ifges(4,4) * t184;
t117 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t184;
t116 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t185;
t103 = Ifges(4,1) * t185 + Ifges(4,5) * qJD(3) + t131;
t102 = Ifges(4,6) * qJD(3) + qJD(2) * t162;
t101 = t104 * qJD(3);
t100 = t105 * qJD(3);
t91 = qJDD(3) * mrSges(4,1) - mrSges(4,3) * t111;
t90 = -qJDD(3) * mrSges(4,2) + mrSges(4,3) * t110;
t89 = Ifges(5,4) * t98;
t76 = -t104 * pkin(4) - t129;
t75 = qJD(3) * mrSges(5,1) + mrSges(5,3) * t99;
t74 = -qJD(3) * mrSges(5,2) + mrSges(5,3) * t98;
t71 = pkin(4) * t100 + t180;
t70 = pkin(3) * t185 - pkin(4) * t99;
t59 = -mrSges(5,1) * t98 - mrSges(5,2) * t99;
t58 = qJDD(3) * mrSges(5,1) - mrSges(5,3) * t64;
t57 = -qJDD(3) * mrSges(5,2) + mrSges(5,3) * t63;
t52 = -t99 * Ifges(5,1) + Ifges(5,5) * qJD(3) + t89;
t51 = t98 * Ifges(5,2) + Ifges(5,6) * qJD(3) - t208;
t40 = mrSges(6,1) * t144 - mrSges(6,3) * t56;
t39 = -mrSges(6,2) * t144 + mrSges(6,3) * t169;
t34 = -pkin(4) * t63 + t88;
t33 = -pkin(7) * t100 + t47;
t31 = -pkin(7) * t101 + t45;
t27 = -qJD(5) * t61 - t100 * t151 - t101 * t149;
t26 = qJD(5) * t60 - t100 * t149 + t101 * t151;
t25 = -mrSges(6,1) * t169 + mrSges(6,2) * t56;
t13 = -mrSges(6,2) * t141 + mrSges(6,3) * t18;
t12 = mrSges(6,1) * t141 - mrSges(6,3) * t17;
t5 = -qJD(5) * t22 - t149 * t33 + t151 * t31;
t4 = qJD(5) * t21 + t149 * t31 + t151 * t33;
t1 = [-t100 * t75 + t101 * t74 + t104 * t58 + t105 * t57 + t60 * t12 + t61 * t13 + t150 * t90 + t152 * t91 + t26 * t39 + t27 * t40 + t214 * qJDD(1) + (-t116 * t150 + t117 * t152) * qJD(3) + m(4) * (t68 * t150 + t152 * t69 + (-t113 * t150 + t114 * t152) * qJD(3)) + m(5) * (-t100 * t41 + t101 * t42 + t104 * t19 + t105 * t20) + m(6) * (t2 * t61 + t26 * t8 + t27 * t7 + t3 * t60) + (-m(4) - m(5) - m(6) - t214) * g(3); (m(4) * t190 + mrSges(4,1) * t110 - mrSges(4,2) * t111) * pkin(2) + t76 * t173 - t129 * t174 + (-t88 * mrSges(5,1) + t20 * mrSges(5,3) + Ifges(5,4) * t64 + t63 * Ifges(5,2) + Ifges(5,6) * t223) * t104 + (t88 * mrSges(5,2) - t19 * mrSges(5,3) + t64 * Ifges(5,1) + Ifges(5,4) * t63 + Ifges(5,5) * t223) * t105 + (-t113 * t188 - t114 * t189 + t224) * mrSges(4,3) + (t152 * t90 + m(4) * (-qJD(3) * t160 + t224) - t116 * t188 - t150 * t91) * pkin(6) + (t133 * t222 + t135 * t221) * g(1) + (t133 * t221 - t135 * t222) * g(2) + t152 * (Ifges(4,4) * t111 + Ifges(4,2) * t110 + Ifges(4,6) * qJDD(3)) / 0.2e1 + t144 * (Ifges(6,5) * t26 + Ifges(6,6) * t27) / 0.2e1 + t101 * t52 / 0.2e1 - t100 * t51 / 0.2e1 + t111 * (Ifges(4,1) * t150 + t198) / 0.2e1 - t119 * t190 + t103 * t188 / 0.2e1 - t102 * t189 / 0.2e1 - t158 * t182 - t117 * t179 + m(5) * (t115 * t180 - t129 * t88 + t19 * t66 + t20 * t67 + t41 * t45 + t42 * t47) + t47 * t74 + t45 * t75 + t65 * (-mrSges(6,1) * t27 + mrSges(6,2) * t26) + t66 * t58 + t67 * t57 + t71 * t25 + t4 * t39 + t5 * t40 + t26 * t24 / 0.2e1 + t21 * t12 + t22 * t13 + t59 * t180 + (-t7 * t26 + t8 * t27) * mrSges(6,3) + (t152 * (-Ifges(4,2) * t150 + t198) + t157) * t182 / 0.2e1 + t169 * (Ifges(6,4) * t26 + Ifges(6,2) * t27) / 0.2e1 + t115 * (mrSges(5,1) * t100 + mrSges(5,2) * t101) + t98 * (Ifges(5,4) * t101 - Ifges(5,2) * t100) / 0.2e1 + qJD(3) * (Ifges(5,5) * t101 - Ifges(5,6) * t100) / 0.2e1 + (-t100 * t42 - t101 * t41) * mrSges(5,3) + (Ifges(5,1) * t101 - Ifges(5,4) * t100) * t215 + m(6) * (t2 * t22 + t21 * t3 + t34 * t76 + t4 * t8 + t5 * t7 + t65 * t71) + (mrSges(6,2) * t34 - mrSges(6,3) * t3 + Ifges(6,1) * t17 + Ifges(6,4) * t18 + Ifges(6,5) * t141) * t61 + (-mrSges(6,1) * t34 + mrSges(6,3) * t2 + Ifges(6,4) * t17 + Ifges(6,2) * t18 + Ifges(6,6) * t141) * t60 + t110 * t162 / 0.2e1 + qJD(3) ^ 2 * t161 / 0.2e1 + Ifges(3,3) * qJDD(2) + t27 * t229 + (Ifges(4,5) * t150 + Ifges(4,6) * t152) * t191 + (Ifges(4,1) * t111 + Ifges(4,4) * t110 + Ifges(4,5) * qJDD(3)) * t210 + (Ifges(6,1) * t26 + Ifges(6,4) * t27) * t217; (Ifges(5,3) + Ifges(4,3)) * qJDD(3) + t99 * (Ifges(5,1) * t98 + t208) / 0.2e1 - t115 * (-mrSges(5,1) * t99 + mrSges(5,2) * t98) + t114 * t116 - t113 * t117 + Ifges(4,6) * t110 + Ifges(4,5) * t111 - m(5) * (t41 * t44 + t42 * t46) + t92 * t12 + t93 * t13 - t46 * t74 - t44 * t75 + Ifges(5,6) * t63 + Ifges(5,5) * t64 - t68 * mrSges(4,2) + t69 * mrSges(4,1) - t70 * t25 + t19 * mrSges(5,1) - t20 * mrSges(5,2) - (Ifges(5,2) * t99 + t52 + t89) * t98 / 0.2e1 + t56 * t229 + (t161 * t201 + t102 * t210 + (t158 + t194 * t210 - t157 / 0.2e1) * qJD(2) + (-m(5) * t115 - t59) * t205 + t160 * mrSges(4,3) - (t131 + t103) * t152 / 0.2e1) * qJD(2) + (t41 * t98 - t42 * t99) * mrSges(5,3) + (t146 * t57 + t147 * t58) * pkin(3) + t235 + (-g(3) * t170 + t2 * t93 + t233 * t7 + t234 * t8 + t3 * t92 - t65 * t70) * m(6) + t233 * t40 + t234 * t39 + (-m(5) * t204 - t231) * g(3) + t232 * (-m(6) * (-pkin(4) * t134 - t205) + mrSges(5,1) * t134 + t200 + mrSges(5,2) * t136 + (mrSges(4,1) + t220) * t150 + t163) + (Ifges(5,5) * t98 + Ifges(5,6) * t99) * t201 + t51 * t215 + (t146 * t20 + t147 * t19) * t220; -t169 * t39 + t56 * t40 - t98 * t74 - t99 * t75 + t173 + t174 + (-t169 * t8 + t56 * t7 + t165 + t34) * m(6) + (-t41 * t99 - t42 * t98 + t165 + t88) * m(5); -g(3) * t166 + t232 * t163 + t23 * t217 - t7 * t39 + t8 * t40 + t235;];
tau = t1;
