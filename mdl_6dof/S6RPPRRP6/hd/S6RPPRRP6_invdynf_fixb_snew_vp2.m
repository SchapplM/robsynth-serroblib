% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP6
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
% Datum: 2019-05-05 15:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:00:11
% EndTime: 2019-05-05 15:00:12
% DurationCPUTime: 0.58s
% Computational Cost: add. (5116->149), mult. (9389->171), div. (0->0), fcn. (4717->6), ass. (0->71)
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t104 = t73 * g(1) - t75 * g(2);
t77 = qJD(1) ^ 2;
t42 = -qJDD(1) * pkin(1) - t77 * qJ(2) + qJDD(2) - t104;
t117 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) + t42;
t89 = -t75 * g(1) - t73 * g(2);
t116 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t89;
t74 = cos(qJ(4));
t102 = qJD(1) * t74;
t72 = sin(qJ(4));
t53 = (pkin(4) * t72 - pkin(8) * t74) * qJD(1);
t76 = qJD(4) ^ 2;
t78 = qJDD(3) + (-pkin(1) - qJ(3)) * t77 + t116;
t32 = -qJDD(1) * pkin(7) + t78;
t88 = t72 * g(3) + t74 * t32;
t19 = -qJDD(4) * pkin(4) - t76 * pkin(8) + t53 * t102 - t88;
t112 = cos(qJ(5));
t71 = sin(qJ(5));
t51 = t71 * qJD(4) + t112 * t102;
t100 = qJD(1) * qJD(4);
t92 = t72 * t100;
t55 = qJDD(1) * t74 - t92;
t23 = qJD(5) * t51 - t112 * qJDD(4) + t55 * t71;
t50 = -t112 * qJD(4) + t71 * t102;
t24 = -t50 * qJD(5) + t71 * qJDD(4) + t112 * t55;
t103 = qJD(1) * t72;
t59 = qJD(5) + t103;
t37 = -mrSges(6,2) * t59 - mrSges(6,3) * t50;
t38 = mrSges(6,1) * t59 - mrSges(6,3) * t51;
t39 = -mrSges(7,1) * t59 + mrSges(7,2) * t51;
t40 = -mrSges(7,2) * t50 + mrSges(7,3) * t59;
t95 = m(7) * (-0.2e1 * qJD(6) * t51 + (t50 * t59 - t24) * qJ(6) + (t51 * t59 + t23) * pkin(5) + t19) + t23 * mrSges(7,1) + t50 * t40;
t115 = m(6) * t19 + t23 * mrSges(6,1) + (t38 - t39) * t51 + (mrSges(6,2) - mrSges(7,3)) * t24 + t50 * t37 + t95;
t114 = -m(3) - m(4);
t27 = pkin(5) * t50 - qJ(6) * t51;
t91 = t74 * t100;
t54 = -qJDD(1) * t72 - t91;
t49 = qJDD(5) - t54;
t58 = t59 ^ 2;
t81 = -t77 * pkin(7) - t117;
t17 = (-t55 + t92) * pkin(8) + (-t54 + t91) * pkin(4) + t81;
t93 = -g(3) * t74 + t72 * t32;
t20 = -pkin(4) * t76 + qJDD(4) * pkin(8) - t53 * t103 + t93;
t84 = t112 * t17 - t71 * t20;
t113 = m(7) * (-t49 * pkin(5) - t58 * qJ(6) + t51 * t27 + qJDD(6) - t84);
t111 = mrSges(3,2) - mrSges(4,3);
t110 = -mrSges(4,2) - mrSges(3,3);
t108 = -mrSges(6,3) - mrSges(7,2);
t107 = t112 * t20 + t71 * t17;
t28 = mrSges(7,1) * t50 - mrSges(7,3) * t51;
t106 = -mrSges(6,1) * t50 - mrSges(6,2) * t51 - t28;
t101 = -m(2) + t114;
t97 = mrSges(2,1) - t111;
t96 = m(7) * (-pkin(5) * t58 + qJ(6) * t49 + 0.2e1 * qJD(6) * t59 - t27 * t50 + t107) + t59 * t39 + t49 * mrSges(7,3);
t11 = m(6) * t84 - t113 + (t37 + t40) * t59 + t106 * t51 + (mrSges(6,1) + mrSges(7,1)) * t49 + t108 * t24;
t52 = (mrSges(5,1) * t72 + mrSges(5,2) * t74) * qJD(1);
t57 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t102;
t9 = m(6) * t107 - t49 * mrSges(6,2) + t106 * t50 + t108 * t23 - t59 * t38 + t96;
t5 = m(5) * t93 - qJDD(4) * mrSges(5,2) + t54 * mrSges(5,3) - qJD(4) * t57 - t52 * t103 - t71 * t11 + t112 * t9;
t56 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t103;
t7 = m(5) * t88 + qJDD(4) * mrSges(5,1) - t55 * mrSges(5,3) + qJD(4) * t56 - t52 * t102 - t115;
t94 = t74 * t5 - t72 * t7;
t90 = m(4) * t78 + qJDD(1) * mrSges(4,2) + t72 * t5 + t74 * t7;
t86 = m(3) * (pkin(1) * t77 - t116) - t90;
t83 = m(5) * t81 - t54 * mrSges(5,1) + t55 * mrSges(5,2) + t57 * t102 + t56 * t103 + t112 * t11 + t71 * t9;
t82 = m(4) * t117 - t83;
t80 = m(3) * t42 + t82;
t2 = m(2) * t104 + (-mrSges(2,2) - t110) * t77 + t97 * qJDD(1) - t80;
t1 = m(2) * t89 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t97 * t77 - t86;
t3 = [-m(1) * g(1) + t1 * t75 - t2 * t73, t1, t114 * g(3) + t94, -m(4) * g(3) + t94, t5, t9, -t23 * mrSges(7,2) - t50 * t28 + t96; -m(1) * g(2) + t1 * t73 + t2 * t75, t2, -qJDD(1) * mrSges(3,3) - t111 * t77 + t86, -t77 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t82, t7, t11, -t24 * mrSges(7,3) - t51 * t39 + t95; (-m(1) + t101) * g(3) + t94, t101 * g(3) + t94, t111 * qJDD(1) + t110 * t77 + t80, -t77 * mrSges(4,3) + t90, t83, t115, -t49 * mrSges(7,1) + t24 * mrSges(7,2) + t51 * t28 - t59 * t40 + t113;];
f_new  = t3;
