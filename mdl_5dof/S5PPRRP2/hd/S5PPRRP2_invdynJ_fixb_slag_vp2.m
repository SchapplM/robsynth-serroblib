% Calculate vector of inverse dynamics joint torques for
% S5PPRRP2
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
%   pkin=[a2,a3,a4,a5,d3,d4,theta1,theta2]';
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
% Datum: 2019-12-05 15:09
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRRP2_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP2_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:08:32
% EndTime: 2019-12-05 15:08:41
% DurationCPUTime: 3.14s
% Computational Cost: add. (1011->252), mult. (2294->327), div. (0->0), fcn. (1489->10), ass. (0->118)
t179 = mrSges(5,1) + mrSges(6,1);
t177 = Ifges(6,4) + Ifges(5,5);
t176 = Ifges(6,6) - Ifges(5,6);
t76 = sin(qJ(4));
t78 = cos(qJ(4));
t100 = mrSges(5,1) * t78 - mrSges(5,2) * t76;
t193 = -mrSges(4,1) - t100;
t128 = t78 * qJD(3);
t112 = mrSges(5,3) * t128;
t113 = mrSges(6,2) * t128;
t49 = qJD(4) * mrSges(6,3) + t113;
t137 = -qJD(4) * mrSges(5,2) + t112 + t49;
t129 = t76 * qJD(3);
t114 = mrSges(6,2) * t129;
t119 = mrSges(5,3) * t129;
t174 = qJD(4) * t179 - t114 - t119;
t82 = t137 * t78 - t174 * t76;
t192 = -qJD(3) * mrSges(4,2) + t82;
t98 = t78 * mrSges(6,1) + t76 * mrSges(6,3);
t191 = -t98 + t193;
t71 = pkin(8) + qJ(3);
t67 = sin(t71);
t160 = g(3) * t67;
t72 = sin(pkin(8));
t74 = cos(pkin(8));
t77 = sin(qJ(3));
t79 = cos(qJ(3));
t183 = -t72 * t77 + t79 * t74;
t29 = t183 * qJD(1);
t92 = pkin(4) * t78 + qJ(5) * t76;
t45 = -pkin(3) - t92;
t13 = qJD(3) * t45 - t29;
t18 = -qJD(3) * pkin(3) - t29;
t97 = t76 * mrSges(6,1) - t78 * mrSges(6,3);
t99 = mrSges(5,1) * t76 + mrSges(5,2) * t78;
t190 = -t13 * t97 - t18 * t99;
t180 = -m(6) - m(5);
t68 = cos(t71);
t188 = t180 * t68;
t131 = qJD(4) * t78;
t125 = qJD(1) * qJD(3);
t184 = qJDD(1) * t79 - t77 * t125;
t185 = qJDD(1) * t77 + t79 * t125;
t9 = t184 * t72 + t185 * t74;
t7 = qJDD(3) * pkin(6) + t9;
t122 = qJD(2) * t131 + t76 * qJDD(2) + t78 * t7;
t132 = qJD(4) * t76;
t88 = t72 * t79 + t77 * t74;
t30 = t88 * qJD(1);
t19 = qJD(3) * pkin(6) + t30;
t3 = -t132 * t19 + t122;
t15 = qJD(2) * t76 + t19 * t78;
t4 = -t15 * qJD(4) + qJDD(2) * t78 - t7 * t76;
t102 = t3 * t78 - t4 * t76;
t151 = t19 * t76;
t14 = qJD(2) * t78 - t151;
t187 = -t14 * t131 - t15 * t132 + t102;
t1 = qJDD(4) * qJ(5) + (qJD(5) - t151) * qJD(4) + t122;
t2 = -qJDD(4) * pkin(4) + qJDD(5) - t4;
t103 = t1 * t78 + t2 * t76;
t11 = -qJD(4) * pkin(4) + qJD(5) - t14;
t12 = qJD(4) * qJ(5) + t15;
t130 = t12 * qJD(4);
t186 = t11 * t131 - t76 * t130 + t103;
t124 = qJD(3) * qJD(4);
t44 = qJDD(3) * t76 + t124 * t78;
t26 = -qJDD(4) * mrSges(6,1) + t44 * mrSges(6,2);
t139 = qJDD(4) * mrSges(5,1) - mrSges(5,3) * t44 - t26;
t43 = -t78 * qJDD(3) + t124 * t76;
t27 = -mrSges(6,2) * t43 + qJDD(4) * mrSges(6,3);
t140 = -qJDD(4) * mrSges(5,2) - mrSges(5,3) * t43 + t27;
t181 = -t139 * t76 + t140 * t78;
t178 = mrSges(5,2) - mrSges(6,3);
t66 = Ifges(5,4) * t128;
t152 = Ifges(6,5) * t78;
t96 = t76 * Ifges(6,1) - t152;
t175 = Ifges(5,1) * t129 + qJD(3) * t96 + qJD(4) * t177 + t66;
t173 = t176 * t76 + t177 * t78;
t40 = t98 * qJD(3);
t169 = qJD(3) * t193 - t40;
t165 = t76 / 0.2e1;
t164 = m(3) + m(4);
t155 = Ifges(5,4) * t76;
t154 = Ifges(5,4) * t78;
t153 = Ifges(6,5) * t76;
t31 = t183 * qJD(3);
t148 = t31 * t76;
t147 = t31 * t78;
t73 = sin(pkin(7));
t145 = t73 * t76;
t144 = t73 * t78;
t75 = cos(pkin(7));
t143 = t75 * t76;
t142 = t75 * t78;
t111 = t164 - t180;
t105 = -t124 / 0.2e1;
t95 = t78 * Ifges(5,2) + t155;
t91 = pkin(4) * t76 - qJ(5) * t78;
t90 = t11 * t76 + t12 * t78;
t89 = -t14 * t76 + t15 * t78;
t84 = t76 * (Ifges(5,1) * t78 - t155);
t83 = t78 * (Ifges(6,3) * t76 + t152);
t10 = t184 * t74 - t185 * t72;
t8 = -qJDD(3) * pkin(3) - t10;
t65 = Ifges(6,5) * t129;
t42 = t91 * qJD(3);
t34 = Ifges(5,6) * qJD(4) + qJD(3) * t95;
t33 = Ifges(6,6) * qJD(4) - Ifges(6,3) * t128 + t65;
t32 = t88 * qJD(3);
t28 = qJD(4) * t91 - qJD(5) * t76;
t23 = t142 * t68 + t145;
t22 = t143 * t68 - t144;
t21 = t144 * t68 - t143;
t20 = t145 * t68 + t142;
t17 = mrSges(5,1) * t43 + mrSges(5,2) * t44;
t16 = mrSges(6,1) * t43 - mrSges(6,3) * t44;
t5 = t43 * pkin(4) - t44 * qJ(5) - qJD(5) * t129 + t8;
t6 = [-(-qJDD(3) * mrSges(4,1) + t16 + t17) * t183 + t169 * t32 + t192 * t31 + (-m(2) - t111) * g(3) + m(6) * (t11 * t148 + t12 * t147 + t13 * t32 - t183 * t5) + m(5) * (-t14 * t148 + t15 * t147 + t18 * t32 - t183 * t8) + m(4) * (t10 * t183 - t29 * t32 + t30 * t31) + (-qJDD(3) * mrSges(4,2) + (-t137 * t76 - t174 * t78) * qJD(4) + m(6) * t186 + m(5) * t187 + m(4) * t9 + t181) * t88 + (m(2) + m(3) * (t72 ^ 2 + t74 ^ 2)) * qJDD(1); t139 * t78 + t140 * t76 + t164 * qJDD(2) + t82 * qJD(4) + m(5) * (qJD(4) * t89 + t3 * t76 + t4 * t78) + m(6) * (qJD(4) * t90 + t1 * t76 - t2 * t78) + (-g(1) * t73 + g(2) * t75) * t111; (g(1) * t75 + g(2) * t73) * (pkin(6) * t188 + (mrSges(4,2) - mrSges(6,2) - mrSges(5,3)) * t68 + (m(5) * pkin(3) - m(6) * t45 - t191) * t67) + (mrSges(4,2) * t67 + (-m(6) * t92 + t191) * t68) * g(3) + (-t160 + t186) * mrSges(6,2) + (-t160 + t187) * mrSges(5,3) + (-m(5) * t8 + g(3) * t188 - t17) * pkin(3) + (-t95 / 0.2e1 + t153 / 0.2e1 + (-Ifges(6,3) - Ifges(5,2) / 0.2e1) * t78 + (Ifges(6,5) - Ifges(5,4)) * t165) * t43 + (-t176 * t78 + t177 * t76) * qJDD(4) / 0.2e1 + ((Ifges(5,1) + Ifges(6,1)) * t44 + t177 * qJDD(4)) * t165 - t8 * t100 - t5 * t98 + t175 * t131 / 0.2e1 + (-m(5) * t18 - m(6) * t13 - t169) * t30 + (t173 * qJD(4) / 0.2e1 - t190) * qJD(4) - t78 * (Ifges(6,5) * t44 + Ifges(6,6) * qJDD(4)) / 0.2e1 + t78 * (Ifges(5,4) * t44 + Ifges(5,6) * qJDD(4)) / 0.2e1 + (-t34 / 0.2e1 + t33 / 0.2e1) * t132 - t28 * t40 + t45 * t16 - t9 * mrSges(4,2) + t10 * mrSges(4,1) + (t180 * t160 + m(5) * ((-t14 * t78 - t15 * t76) * qJD(4) + t102) + m(6) * ((t11 * t78 - t12 * t76) * qJD(4) + t103) - t137 * t132 - t174 * t131 + t181) * pkin(6) + (t76 * (Ifges(6,1) * t78 + t153) + t78 * (-Ifges(5,2) * t76 + t154) + t84) * t124 / 0.2e1 + (t76 * Ifges(5,1) + t154 + t96) * t44 / 0.2e1 + t83 * t105 + (-m(5) * t89 - m(6) * t90 - t192) * t29 + m(6) * (t13 * t28 + t45 * t5) + Ifges(4,3) * qJDD(3); t34 * t129 / 0.2e1 - t11 * t113 + (Ifges(6,2) + Ifges(5,3)) * qJDD(4) + t176 * t43 + t177 * t44 + (t178 * t23 + t179 * t22) * g(1) + (t178 * t21 + t179 * t20) * g(2) + t173 * t105 + (-m(6) * t12 + t112 - t137) * t14 + (-m(6) * t11 + t119 + t174) * t15 - (-Ifges(5,2) * t129 + t175 + t66) * t128 / 0.2e1 + (t97 + t99) * t160 + qJD(5) * t49 + t42 * t40 + qJ(5) * t27 - pkin(4) * t26 + t12 * t114 + t1 * mrSges(6,3) - t2 * mrSges(6,1) - t3 * mrSges(5,2) + t4 * mrSges(5,1) - (Ifges(6,1) * t128 + t33 + t65) * t129 / 0.2e1 + (-pkin(4) * t2 + qJ(5) * t1 + qJD(5) * t12 - g(2) * (-pkin(4) * t20 + qJ(5) * t21) - g(1) * (-pkin(4) * t22 + qJ(5) * t23) + t91 * t160 - t13 * t42) * m(6) + ((-t84 / 0.2e1 + t83 / 0.2e1) * qJD(3) + t190) * qJD(3); -t40 * t129 - qJD(4) * t49 + (-g(1) * t22 - g(2) * t20 + t129 * t13 - t160 * t76 - t130 + t2) * m(6) + t26;];
tau = t6;
