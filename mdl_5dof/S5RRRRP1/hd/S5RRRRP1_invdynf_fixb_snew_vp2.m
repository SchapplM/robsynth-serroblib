% Calculate vector of cutting forces with Newton-Euler
% S5RRRRP1
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4]';
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
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:46
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:44:49
% EndTime: 2019-12-05 18:44:53
% DurationCPUTime: 1.56s
% Computational Cost: add. (16117->162), mult. (35361->208), div. (0->0), fcn. (24282->8), ass. (0->77)
t101 = qJD(1) * qJD(2);
t79 = sin(qJ(2));
t83 = cos(qJ(2));
t65 = t79 * qJDD(1) + t83 * t101;
t66 = t83 * qJDD(1) - t79 * t101;
t103 = qJD(1) * t79;
t67 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t103;
t102 = qJD(1) * t83;
t68 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t102;
t85 = qJD(1) ^ 2;
t78 = sin(qJ(3));
t82 = cos(qJ(3));
t60 = (t78 * t83 + t79 * t82) * qJD(1);
t40 = -t60 * qJD(3) - t78 * t65 + t82 * t66;
t59 = (-t78 * t79 + t82 * t83) * qJD(1);
t41 = t59 * qJD(3) + t82 * t65 + t78 * t66;
t75 = qJD(2) + qJD(3);
t54 = -t75 * mrSges(4,2) + t59 * mrSges(4,3);
t55 = t75 * mrSges(4,1) - t60 * mrSges(4,3);
t69 = qJD(2) * pkin(2) - pkin(7) * t103;
t76 = t83 ^ 2;
t80 = sin(qJ(1));
t84 = cos(qJ(1));
t97 = t80 * g(1) - t84 * g(2);
t90 = -qJDD(1) * pkin(1) - t97;
t88 = -t66 * pkin(2) + t69 * t103 + (-pkin(7) * t76 - pkin(6)) * t85 + t90;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t52 = t77 * t59 + t81 * t60;
t25 = -t52 * qJD(4) + t81 * t40 - t77 * t41;
t51 = t81 * t59 - t77 * t60;
t26 = t51 * qJD(4) + t77 * t40 + t81 * t41;
t72 = qJD(4) + t75;
t43 = -t72 * mrSges(6,2) + t51 * mrSges(6,3);
t44 = -t72 * mrSges(5,2) + t51 * mrSges(5,3);
t47 = t72 * mrSges(5,1) - t52 * mrSges(5,3);
t56 = t75 * pkin(3) - t60 * pkin(8);
t58 = t59 ^ 2;
t87 = -t40 * pkin(3) - t58 * pkin(8) + t60 * t56 + t88;
t45 = t72 * pkin(4) - t52 * qJ(5);
t46 = t72 * mrSges(6,1) - t52 * mrSges(6,3);
t50 = t51 ^ 2;
t98 = m(6) * (-t25 * pkin(4) - t50 * qJ(5) + t52 * t45 + qJDD(5) + t87) + t26 * mrSges(6,2) + t52 * t46;
t89 = m(5) * t87 + t26 * mrSges(5,2) - (t44 + t43) * t51 - (mrSges(5,1) + mrSges(6,1)) * t25 + t52 * t47 + t98;
t86 = m(4) * t88 - t40 * mrSges(4,1) + t41 * mrSges(4,2) - t59 * t54 + t60 * t55 + t89;
t117 = (t79 * t67 - t83 * t68) * qJD(1) - t66 * mrSges(3,1) + t65 * mrSges(3,2) + m(3) * (-t85 * pkin(6) + t90) + t86;
t93 = -t84 * g(1) - t80 * g(2);
t62 = -t85 * pkin(1) + qJDD(1) * pkin(6) + t93;
t108 = t79 * t62;
t74 = qJDD(2) + qJDD(3);
t113 = pkin(2) * t85;
t35 = qJDD(2) * pkin(2) - t65 * pkin(7) - t108 + (pkin(7) * t101 + t79 * t113 - g(3)) * t83;
t96 = -t79 * g(3) + t83 * t62;
t36 = t66 * pkin(7) - qJD(2) * t69 - t76 * t113 + t96;
t94 = t82 * t35 - t78 * t36;
t17 = (t59 * t75 - t41) * pkin(8) + (t59 * t60 + t74) * pkin(3) + t94;
t105 = t78 * t35 + t82 * t36;
t19 = -t58 * pkin(3) + t40 * pkin(8) - t75 * t56 + t105;
t106 = t77 * t17 + t81 * t19;
t31 = -t51 * mrSges(5,1) + t52 * mrSges(5,2);
t71 = qJDD(4) + t74;
t30 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t99 = m(6) * (-t50 * pkin(4) + t25 * qJ(5) + 0.2e1 * qJD(5) * t51 - t72 * t45 + t106) + t25 * mrSges(6,3) + t51 * t30;
t10 = m(5) * t106 + t25 * mrSges(5,3) + t51 * t31 + (-t47 - t46) * t72 + (-mrSges(5,2) - mrSges(6,2)) * t71 + t99;
t53 = -t59 * mrSges(4,1) + t60 * mrSges(4,2);
t95 = t81 * t17 - t77 * t19;
t100 = m(6) * (-0.2e1 * qJD(5) * t52 + (t51 * t72 - t26) * qJ(5) + (t51 * t52 + t71) * pkin(4) + t95) + t72 * t43 + t71 * mrSges(6,1);
t8 = m(5) * t95 + t71 * mrSges(5,1) + t72 * t44 + (-t31 - t30) * t52 + (-mrSges(5,3) - mrSges(6,3)) * t26 + t100;
t6 = m(4) * t94 + t74 * mrSges(4,1) - t41 * mrSges(4,3) + t77 * t10 - t60 * t53 + t75 * t54 + t81 * t8;
t64 = (-mrSges(3,1) * t83 + mrSges(3,2) * t79) * qJD(1);
t7 = m(4) * t105 - t74 * mrSges(4,2) + t40 * mrSges(4,3) + t81 * t10 + t59 * t53 - t75 * t55 - t77 * t8;
t4 = m(3) * (-t83 * g(3) - t108) - t65 * mrSges(3,3) + qJDD(2) * mrSges(3,1) - t64 * t103 + qJD(2) * t68 + t78 * t7 + t82 * t6;
t5 = m(3) * t96 - qJDD(2) * mrSges(3,2) + t66 * mrSges(3,3) - qJD(2) * t67 + t64 * t102 - t78 * t6 + t82 * t7;
t115 = t83 * t4 + t79 * t5;
t9 = m(2) * t97 + qJDD(1) * mrSges(2,1) - t85 * mrSges(2,2) - t117;
t1 = m(2) * t93 - t85 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t79 * t4 + t83 * t5;
t2 = [-m(1) * g(1) + t84 * t1 - t80 * t9, t1, t5, t7, t10, -t71 * mrSges(6,2) - t72 * t46 + t99; -m(1) * g(2) + t80 * t1 + t84 * t9, t9, t4, t6, t8, -t26 * mrSges(6,3) - t52 * t30 + t100; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t117, t86, t89, -t25 * mrSges(6,1) - t51 * t43 + t98;];
f_new = t2;
