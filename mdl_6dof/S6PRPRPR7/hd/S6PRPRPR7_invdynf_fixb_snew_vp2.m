% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR7
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-05-04 23:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:18:56
% EndTime: 2019-05-04 23:18:58
% DurationCPUTime: 0.91s
% Computational Cost: add. (7632->150), mult. (14304->188), div. (0->0), fcn. (8096->10), ass. (0->82)
t85 = qJD(2) ^ 2;
t120 = -pkin(2) - pkin(8);
t105 = t120 * t85;
t82 = cos(qJ(4));
t109 = t82 * qJD(2);
t79 = sin(qJ(4));
t110 = qJD(2) * t79;
t114 = mrSges(5,1) - mrSges(6,2);
t108 = qJD(2) * qJD(4);
t103 = t82 * t108;
t54 = t79 * qJDD(2) + t103;
t68 = t79 * t108;
t55 = t82 * qJDD(2) - t68;
t58 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t110;
t59 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t109;
t74 = sin(pkin(10));
t76 = cos(pkin(10));
t56 = t74 * g(1) - t76 * g(2);
t77 = cos(pkin(6));
t117 = t56 * t77;
t57 = -t76 * g(1) - t74 * g(2);
t73 = -g(3) + qJDD(1);
t75 = sin(pkin(6));
t80 = sin(qJ(2));
t83 = cos(qJ(2));
t123 = -t80 * t57 + (t73 * t75 + t117) * t83;
t88 = -t85 * qJ(3) + qJDD(3) - t123;
t26 = t120 * qJDD(2) + t88;
t39 = -t75 * t56 + t77 * t73;
t102 = t82 * t26 - t79 * t39;
t51 = (pkin(4) * t79 - qJ(5) * t82) * qJD(2);
t84 = qJD(4) ^ 2;
t20 = -qJDD(4) * pkin(4) - t84 * qJ(5) + t51 * t109 + qJDD(5) - t102;
t17 = (t79 * t82 * t85 - qJDD(4)) * pkin(9) + (t55 + t68) * pkin(5) + t20;
t63 = pkin(5) * t109 - qJD(4) * pkin(9);
t72 = t79 ^ 2;
t121 = -2 * qJD(5);
t116 = t75 * t80;
t106 = t73 * t116 + t80 * t117 + t83 * t57;
t98 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t106;
t90 = pkin(4) * t103 + t109 * t121 + t98 + (-t55 + t68) * qJ(5);
t18 = -t63 * t109 + (pkin(4) + pkin(9)) * t54 + (-pkin(5) * t72 + t120) * t85 + t90;
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t49 = -t78 * qJD(4) + t81 * t110;
t32 = t49 * qJD(6) + t81 * qJDD(4) + t78 * t54;
t50 = t81 * qJD(4) + t78 * t110;
t33 = -t49 * mrSges(7,1) + t50 * mrSges(7,2);
t66 = qJD(6) + t109;
t36 = -t66 * mrSges(7,2) + t49 * mrSges(7,3);
t46 = qJDD(6) + t55;
t13 = m(7) * (t81 * t17 - t78 * t18) - t32 * mrSges(7,3) + t46 * mrSges(7,1) - t50 * t33 + t66 * t36;
t31 = -t50 * qJD(6) - t78 * qJDD(4) + t81 * t54;
t37 = t66 * mrSges(7,1) - t50 * mrSges(7,3);
t14 = m(7) * (t78 * t17 + t81 * t18) + t31 * mrSges(7,3) - t46 * mrSges(7,2) + t49 * t33 - t66 * t37;
t60 = mrSges(6,1) * t110 - qJD(4) * mrSges(6,3);
t61 = mrSges(6,1) * t109 + qJD(4) * mrSges(6,2);
t87 = -(t79 * t60 + t82 * t61) * qJD(2) - t78 * t13 + t81 * t14 + m(6) * (t54 * pkin(4) + t105 + t90) - t55 * mrSges(6,3);
t86 = t114 * t54 + m(5) * (t105 + t98) + t58 * t110 + t59 * t109 + t55 * mrSges(5,2) + t87;
t124 = m(4) * (t85 * pkin(2) - t98) - t86;
t113 = -mrSges(3,2) + mrSges(4,3);
t115 = mrSges(3,1) - mrSges(4,2);
t52 = (-mrSges(6,2) * t79 - mrSges(6,3) * t82) * qJD(2);
t101 = qJD(2) * (-t52 - (mrSges(5,1) * t79 + mrSges(5,2) * t82) * qJD(2));
t111 = t79 * t26 + t82 * t39;
t112 = -mrSges(5,3) - mrSges(6,1);
t89 = -t84 * pkin(4) + qJDD(4) * qJ(5) - t51 * t110 + t111;
t92 = -t31 * mrSges(7,1) - t49 * t36 + m(7) * (-t72 * t85 * pkin(9) - t54 * pkin(5) + ((2 * qJD(5)) + t63) * qJD(4) + t89) + t32 * mrSges(7,2) + t50 * t37;
t91 = -m(6) * (qJD(4) * t121 - t89) + t92;
t11 = m(5) * t111 + t112 * t54 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) + (-t59 + t61) * qJD(4) + t79 * t101 + t91;
t93 = -m(6) * t20 - t81 * t13 - t78 * t14;
t9 = m(5) * t102 + t112 * t55 + t114 * qJDD(4) + (t58 - t60) * qJD(4) + t82 * t101 + t93;
t94 = -m(4) * (-qJDD(2) * pkin(2) + t88) - t79 * t11 - t82 * t9;
t4 = m(3) * t123 + t115 * qJDD(2) + t113 * t85 + t94;
t118 = t4 * t83;
t99 = m(4) * t39 + t82 * t11 - t79 * t9;
t6 = m(3) * t39 + t99;
t8 = m(3) * t106 + t113 * qJDD(2) - t115 * t85 - t124;
t104 = m(2) * t73 + t8 * t116 + t75 * t118 + t77 * t6;
t2 = m(2) * t57 - t80 * t4 + t83 * t8;
t1 = m(2) * t56 - t75 * t6 + (t8 * t80 + t118) * t77;
t3 = [-m(1) * g(1) - t74 * t1 + t76 * t2, t2, t8, t99, t11, -t54 * mrSges(6,2) + t87, t14; -m(1) * g(2) + t76 * t1 + t74 * t2, t1, t4, -t85 * mrSges(4,2) - qJDD(2) * mrSges(4,3) + t124, t9, t54 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t61 + t52 * t110 - t91, t13; -m(1) * g(3) + t104, t104, t6, qJDD(2) * mrSges(4,2) - t85 * mrSges(4,3) - t94, t86, t55 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t60 + t52 * t109 - t93, t92;];
f_new  = t3;
