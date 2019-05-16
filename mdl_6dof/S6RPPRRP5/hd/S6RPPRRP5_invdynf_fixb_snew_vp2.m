% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP5
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
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5]';
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
% Datum: 2019-05-05 14:59
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:57:26
% EndTime: 2019-05-05 14:57:28
% DurationCPUTime: 0.60s
% Computational Cost: add. (5182->149), mult. (9605->171), div. (0->0), fcn. (4891->6), ass. (0->70)
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t107 = t74 * g(1) - t77 * g(2);
t79 = qJD(1) ^ 2;
t45 = -qJDD(1) * pkin(1) - t79 * qJ(2) + qJDD(2) - t107;
t116 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t45;
t90 = -t77 * g(1) - t74 * g(2);
t115 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t90;
t76 = cos(qJ(4));
t105 = qJD(1) * t76;
t73 = sin(qJ(4));
t57 = (pkin(4) * t73 - pkin(8) * t76) * qJD(1);
t78 = qJD(4) ^ 2;
t80 = qJDD(3) + (-pkin(1) - qJ(3)) * t79 + t115;
t34 = -qJDD(1) * pkin(7) + t80;
t89 = g(3) * t73 + t34 * t76;
t20 = -qJDD(4) * pkin(4) - pkin(8) * t78 + t57 * t105 - t89;
t72 = sin(qJ(5));
t75 = cos(qJ(5));
t55 = qJD(4) * t72 + t75 * t105;
t103 = qJD(1) * qJD(4);
t94 = t73 * t103;
t59 = qJDD(1) * t76 - t94;
t26 = -qJD(5) * t55 + qJDD(4) * t75 - t59 * t72;
t54 = qJD(4) * t75 - t72 * t105;
t27 = qJD(5) * t54 + qJDD(4) * t72 + t59 * t75;
t106 = qJD(1) * t73;
t62 = qJD(5) + t106;
t39 = -mrSges(7,2) * t62 + mrSges(7,3) * t54;
t40 = -mrSges(6,2) * t62 + mrSges(6,3) * t54;
t43 = mrSges(6,1) * t62 - mrSges(6,3) * t55;
t41 = pkin(5) * t62 - qJ(6) * t55;
t42 = mrSges(7,1) * t62 - mrSges(7,3) * t55;
t52 = t54 ^ 2;
t97 = m(7) * (-pkin(5) * t26 - qJ(6) * t52 + t41 * t55 + qJDD(6) + t20) + t27 * mrSges(7,2) + t55 * t42;
t114 = m(6) * t20 + t27 * mrSges(6,2) - (t40 + t39) * t54 - (mrSges(6,1) + mrSges(7,1)) * t26 + t55 * t43 + t97;
t113 = -m(3) - m(4);
t111 = mrSges(3,2) - mrSges(4,3);
t110 = -mrSges(4,2) - mrSges(3,3);
t93 = t76 * t103;
t58 = -qJDD(1) * t73 - t93;
t83 = -pkin(7) * t79 - t116;
t18 = (-t59 + t94) * pkin(8) + (-t58 + t93) * pkin(4) + t83;
t96 = -g(3) * t76 + t73 * t34;
t21 = -pkin(4) * t78 + qJDD(4) * pkin(8) - t57 * t106 + t96;
t109 = t72 * t18 + t75 * t21;
t104 = -m(2) + t113;
t100 = mrSges(2,1) - t111;
t53 = qJDD(5) - t58;
t92 = t75 * t18 - t72 * t21;
t99 = m(7) * (-0.2e1 * qJD(6) * t55 + (t54 * t62 - t27) * qJ(6) + (t54 * t55 + t53) * pkin(5) + t92) + t62 * t39 + t53 * mrSges(7,1);
t30 = -mrSges(7,1) * t54 + mrSges(7,2) * t55;
t98 = m(7) * (-pkin(5) * t52 + qJ(6) * t26 + 0.2e1 * qJD(6) * t54 - t41 * t62 + t109) + t54 * t30 + t26 * mrSges(7,3);
t56 = (mrSges(5,1) * t73 + mrSges(5,2) * t76) * qJD(1);
t60 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t106;
t11 = m(5) * t89 + qJDD(4) * mrSges(5,1) - t59 * mrSges(5,3) + qJD(4) * t60 - t56 * t105 - t114;
t61 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t105;
t31 = -mrSges(6,1) * t54 + mrSges(6,2) * t55;
t7 = m(6) * t92 + t53 * mrSges(6,1) + t62 * t40 + (-t31 - t30) * t55 + (-mrSges(6,3) - mrSges(7,3)) * t27 + t99;
t9 = m(6) * t109 + t26 * mrSges(6,3) + t54 * t31 + (-t43 - t42) * t62 + (-mrSges(6,2) - mrSges(7,2)) * t53 + t98;
t5 = m(5) * t96 - qJDD(4) * mrSges(5,2) + t58 * mrSges(5,3) - qJD(4) * t61 - t56 * t106 - t72 * t7 + t75 * t9;
t95 = -t73 * t11 + t76 * t5;
t91 = m(4) * t80 + qJDD(1) * mrSges(4,2) + t76 * t11 + t73 * t5;
t87 = m(3) * (pkin(1) * t79 - t115) - t91;
t85 = m(5) * t83 - t58 * mrSges(5,1) + t59 * mrSges(5,2) + t61 * t105 + t60 * t106 + t75 * t7 + t72 * t9;
t84 = m(4) * t116 - t85;
t82 = m(3) * t45 + t84;
t2 = m(2) * t107 + (-mrSges(2,2) - t110) * t79 + t100 * qJDD(1) - t82;
t1 = m(2) * t90 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t100 * t79 - t87;
t3 = [-m(1) * g(1) + t1 * t77 - t2 * t74, t1, t113 * g(3) + t95, -m(4) * g(3) + t95, t5, t9, -t53 * mrSges(7,2) - t62 * t42 + t98; -m(1) * g(2) + t1 * t74 + t2 * t77, t2, -qJDD(1) * mrSges(3,3) - t111 * t79 + t87, -t79 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t84, t11, t7, -t27 * mrSges(7,3) - t55 * t30 + t99; (-m(1) + t104) * g(3) + t95, t104 * g(3) + t95, t111 * qJDD(1) + t110 * t79 + t82, -t79 * mrSges(4,3) + t91, t85, t114, -t26 * mrSges(7,1) - t54 * t39 + t97;];
f_new  = t3;
