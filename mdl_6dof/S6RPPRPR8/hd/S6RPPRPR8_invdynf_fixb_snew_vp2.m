% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR8
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-05 14:39
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:37:34
% EndTime: 2019-05-05 14:37:37
% DurationCPUTime: 1.08s
% Computational Cost: add. (9895->165), mult. (22494->198), div. (0->0), fcn. (14672->8), ass. (0->86)
t121 = -pkin(1) - qJ(3);
t81 = sin(qJ(1));
t84 = cos(qJ(1));
t107 = t81 * g(1) - t84 * g(2);
t86 = qJD(1) ^ 2;
t96 = -t86 * qJ(2) + qJDD(2) - t107;
t133 = -(2 * qJD(1) * qJD(3)) + qJDD(1) * t121 + t96;
t103 = -t84 * g(1) - t81 * g(2);
t131 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t103;
t79 = cos(pkin(9));
t113 = qJDD(1) * t79;
t78 = sin(pkin(9));
t114 = qJDD(1) * t78;
t126 = sin(qJ(4));
t108 = t78 * t126;
t83 = cos(qJ(4));
t63 = (t79 * t83 - t108) * qJD(1);
t115 = t63 * qJD(4);
t99 = t126 * t79 + t78 * t83;
t42 = qJDD(1) * t99 + t115;
t54 = t63 * pkin(5) - qJD(4) * pkin(8);
t62 = t99 * qJD(1);
t61 = t62 ^ 2;
t116 = t62 * qJD(4);
t130 = -2 * qJD(5);
t43 = -qJDD(1) * t108 + t113 * t83 - t116;
t75 = t78 ^ 2;
t117 = t79 ^ 2 + t75;
t94 = qJDD(3) + t131;
t90 = pkin(3) * t114 + (-pkin(7) * t117 + t121) * t86 + t94;
t87 = pkin(4) * t115 + t63 * t130 + (-t43 + t116) * qJ(5) + t90;
t11 = -t61 * pkin(5) - t63 * t54 + (pkin(4) + pkin(8)) * t42 + t87;
t111 = t78 * g(3) + t133 * t79;
t127 = pkin(3) * t86;
t28 = (-pkin(7) * qJDD(1) - t127 * t78) * t79 + t111;
t104 = -t79 * g(3) + t133 * t78;
t29 = -pkin(7) * t114 - t127 * t75 + t104;
t105 = -t126 * t29 + t83 * t28;
t36 = t62 * pkin(4) - t63 * qJ(5);
t85 = qJD(4) ^ 2;
t18 = -qJDD(4) * pkin(4) - t85 * qJ(5) + t63 * t36 + qJDD(5) - t105;
t12 = (t62 * t63 - qJDD(4)) * pkin(8) + (t43 + t116) * pkin(5) + t18;
t80 = sin(qJ(6));
t82 = cos(qJ(6));
t46 = t82 * qJD(4) + t80 * t62;
t21 = -t46 * qJD(6) - t80 * qJDD(4) + t82 * t42;
t45 = -t80 * qJD(4) + t82 * t62;
t24 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t59 = qJD(6) + t63;
t31 = t59 * mrSges(7,1) - t46 * mrSges(7,3);
t41 = qJDD(6) + t43;
t10 = m(7) * (t82 * t11 + t80 * t12) + t21 * mrSges(7,3) - t41 * mrSges(7,2) + t45 * t24 - t59 * t31;
t53 = t63 * mrSges(6,1) + qJD(4) * mrSges(6,2);
t22 = t45 * qJD(6) + t82 * qJDD(4) + t80 * t42;
t30 = -t59 * mrSges(7,2) + t45 * mrSges(7,3);
t9 = m(7) * (-t80 * t11 + t82 * t12) - t22 * mrSges(7,3) + t41 * mrSges(7,1) - t46 * t24 + t59 * t30;
t100 = -t82 * t10 + t80 * t9 - m(6) * (t42 * pkin(4) + t87) + t63 * t53 + t43 * mrSges(6,3);
t52 = t62 * mrSges(6,1) - qJD(4) * mrSges(6,3);
t118 = -qJD(4) * mrSges(5,2) - t62 * mrSges(5,3) - t52;
t124 = mrSges(5,1) - mrSges(6,2);
t51 = qJD(4) * mrSges(5,1) - t63 * mrSges(5,3);
t89 = m(5) * t90 + t43 * mrSges(5,2) + t118 * t62 + t124 * t42 + t63 * t51 - t100;
t88 = m(4) * (t121 * t86 + t94) + mrSges(4,1) * t114 + mrSges(4,2) * t113 + t89;
t132 = m(3) * (t86 * pkin(1) - t131) - t88;
t106 = t117 * mrSges(4,3);
t129 = -m(2) - m(3);
t125 = mrSges(2,1) - mrSges(3,2);
t123 = -mrSges(2,2) + mrSges(3,3);
t122 = -mrSges(5,3) - mrSges(6,1);
t120 = t126 * t28 + t83 * t29;
t38 = -t62 * mrSges(6,2) - t63 * mrSges(6,3);
t119 = -t62 * mrSges(5,1) - t63 * mrSges(5,2) - t38;
t101 = -qJDD(1) * mrSges(4,3) - t86 * (mrSges(4,1) * t78 + mrSges(4,2) * t79);
t97 = -m(6) * t18 - t80 * t10 - t82 * t9;
t6 = m(5) * t105 + qJD(4) * t118 + qJDD(4) * t124 + t119 * t63 + t122 * t43 + t97;
t92 = -t85 * pkin(4) + qJDD(4) * qJ(5) - t62 * t36 + t120;
t95 = -t21 * mrSges(7,1) - t45 * t30 + m(7) * (-t42 * pkin(5) - t61 * pkin(8) + ((2 * qJD(5)) + t54) * qJD(4) + t92) + t22 * mrSges(7,2) + t46 * t31;
t91 = -m(6) * (qJD(4) * t130 - t92) + t95;
t7 = m(5) * t120 + t119 * t62 + t122 * t42 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) + (-t51 + t53) * qJD(4) + t91;
t3 = m(4) * t111 + t101 * t79 + t126 * t7 + t83 * t6;
t4 = m(4) * t104 + t101 * t78 - t126 * t6 + t83 * t7;
t109 = -t78 * t3 + t79 * t4;
t98 = -m(3) * (-qJDD(1) * pkin(1) + t96) - t79 * t3 - t78 * t4;
t5 = m(2) * t103 + (-t106 - t125) * t86 + t123 * qJDD(1) - t132;
t1 = m(2) * t107 + qJDD(1) * t125 + t123 * t86 + t98;
t2 = [-m(1) * g(1) - t81 * t1 + t84 * t5, t5, -m(3) * g(3) + t109, t4, t7, -t42 * mrSges(6,2) - t62 * t52 - t100, t10; -m(1) * g(2) + t84 * t1 + t81 * t5, t1, -qJDD(1) * mrSges(3,3) + (-mrSges(3,2) + t106) * t86 + t132, t3, t6, t42 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t53 + t62 * t38 - t91, t9; (-m(1) + t129) * g(3) + t109, g(3) * t129 + t109, qJDD(1) * mrSges(3,2) - t86 * mrSges(3,3) - t98, -t106 * t86 + t88, t89, t43 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t52 + t63 * t38 - t97, t95;];
f_new  = t2;
