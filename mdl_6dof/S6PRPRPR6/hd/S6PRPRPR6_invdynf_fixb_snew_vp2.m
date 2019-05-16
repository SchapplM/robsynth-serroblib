% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR6
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta5]';
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
% Datum: 2019-05-04 23:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:07:49
% EndTime: 2019-05-04 23:07:53
% DurationCPUTime: 1.60s
% Computational Cost: add. (21017->149), mult. (41677->198), div. (0->0), fcn. (27396->12), ass. (0->87)
t82 = sin(pkin(6));
t88 = sin(qJ(2));
t118 = t82 * t88;
t81 = sin(pkin(10));
t84 = cos(pkin(10));
t70 = t81 * g(1) - t84 * g(2);
t85 = cos(pkin(6));
t119 = t70 * t85;
t71 = -t84 * g(1) - t81 * g(2);
t79 = -g(3) + qJDD(1);
t91 = cos(qJ(2));
t109 = t79 * t118 + t88 * t119 + t91 * t71;
t123 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t109;
t122 = (t79 * t82 + t119) * t91 - t88 * t71;
t121 = -pkin(2) - pkin(8);
t87 = sin(qJ(4));
t113 = t87 * qJD(2);
t93 = qJD(2) ^ 2;
t95 = -t93 * qJ(3) + qJDD(3) - t122;
t36 = t121 * qJDD(2) + t95;
t52 = -t82 * t70 + t85 * t79;
t90 = cos(qJ(4));
t115 = t87 * t36 + t90 * t52;
t66 = (pkin(4) * t87 - qJ(5) * t90) * qJD(2);
t92 = qJD(4) ^ 2;
t24 = -t92 * pkin(4) + qJDD(4) * qJ(5) - t66 * t113 + t115;
t112 = qJD(2) * qJD(4);
t106 = t90 * t112;
t107 = t87 * t112;
t68 = t87 * qJDD(2) + t106;
t69 = t90 * qJDD(2) - t107;
t97 = t121 * t93 - t123;
t27 = (-t69 + t107) * qJ(5) + (t68 + t106) * pkin(4) + t97;
t114 = qJD(2) * t90;
t80 = sin(pkin(11));
t83 = cos(pkin(11));
t61 = t80 * qJD(4) + t83 * t114;
t102 = -0.2e1 * qJD(5) * t61 - t80 * t24 + t83 * t27;
t50 = t80 * qJDD(4) + t83 * t69;
t60 = t83 * qJD(4) - t80 * t114;
t18 = (t60 * t113 - t50) * pkin(9) + (t60 * t61 + t68) * pkin(5) + t102;
t110 = 0.2e1 * qJD(5) * t60 + t83 * t24 + t80 * t27;
t49 = t83 * qJDD(4) - t80 * t69;
t51 = pkin(5) * t113 - t61 * pkin(9);
t59 = t60 ^ 2;
t19 = -t59 * pkin(5) + t49 * pkin(9) - t51 * t113 + t110;
t86 = sin(qJ(6));
t89 = cos(qJ(6));
t41 = t89 * t60 - t86 * t61;
t30 = t41 * qJD(6) + t86 * t49 + t89 * t50;
t42 = t86 * t60 + t89 * t61;
t32 = -t41 * mrSges(7,1) + t42 * mrSges(7,2);
t75 = qJD(6) + t113;
t39 = -t75 * mrSges(7,2) + t41 * mrSges(7,3);
t63 = qJDD(6) + t68;
t16 = m(7) * (t89 * t18 - t86 * t19) - t30 * mrSges(7,3) + t63 * mrSges(7,1) - t42 * t32 + t75 * t39;
t29 = -t42 * qJD(6) + t89 * t49 - t86 * t50;
t40 = t75 * mrSges(7,1) - t42 * mrSges(7,3);
t17 = m(7) * (t86 * t18 + t89 * t19) + t29 * mrSges(7,3) - t63 * mrSges(7,2) + t41 * t32 - t75 * t40;
t43 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t47 = -mrSges(6,2) * t113 + t60 * mrSges(6,3);
t13 = m(6) * t102 + t68 * mrSges(6,1) - t50 * mrSges(6,3) + t47 * t113 + t89 * t16 + t86 * t17 - t61 * t43;
t48 = mrSges(6,1) * t113 - t61 * mrSges(6,3);
t14 = m(6) * t110 - t68 * mrSges(6,2) + t49 * mrSges(6,3) - t48 * t113 - t86 * t16 + t89 * t17 + t60 * t43;
t67 = (mrSges(5,1) * t87 + mrSges(5,2) * t90) * qJD(2);
t73 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t114;
t10 = m(5) * t115 - qJDD(4) * mrSges(5,2) - t68 * mrSges(5,3) - qJD(4) * t73 - t67 * t113 - t80 * t13 + t83 * t14;
t105 = t90 * t36 - t87 * t52;
t72 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t113;
t23 = -qJDD(4) * pkin(4) - t92 * qJ(5) + t66 * t114 + qJDD(5) - t105;
t99 = t29 * mrSges(7,1) + t41 * t39 - m(7) * (-t49 * pkin(5) - t59 * pkin(9) + t61 * t51 + t23) - t30 * mrSges(7,2) - t42 * t40;
t94 = m(6) * t23 - t49 * mrSges(6,1) + t50 * mrSges(6,2) - t60 * t47 + t61 * t48 - t99;
t15 = m(5) * t105 + qJDD(4) * mrSges(5,1) - t69 * mrSges(5,3) + qJD(4) * t72 - t67 * t114 - t94;
t100 = -m(4) * (-qJDD(2) * pkin(2) + t95) - t87 * t10 - t90 * t15;
t116 = (-mrSges(3,2) + mrSges(4,3));
t117 = mrSges(3,1) - mrSges(4,2);
t4 = m(3) * t122 + t117 * qJDD(2) + (t116 * t93) + t100;
t120 = t4 * t91;
t103 = m(4) * t52 + t90 * t10 - t87 * t15;
t6 = m(3) * t52 + t103;
t98 = m(5) * t97 + t68 * mrSges(5,1) + t69 * mrSges(5,2) + t72 * t113 + t73 * t114 + t83 * t13 + t80 * t14;
t96 = -m(4) * (t93 * pkin(2) + t123) + t98;
t8 = m(3) * t109 + t116 * qJDD(2) - t117 * t93 + t96;
t108 = m(2) * t79 + t8 * t118 + t82 * t120 + t85 * t6;
t2 = m(2) * t71 - t88 * t4 + t91 * t8;
t1 = m(2) * t70 - t82 * t6 + (t8 * t88 + t120) * t85;
t3 = [-m(1) * g(1) - t81 * t1 + t84 * t2, t2, t8, t103, t10, t14, t17; -m(1) * g(2) + t84 * t1 + t81 * t2, t1, t4, -(t93 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t96, t15, t13, t16; -m(1) * g(3) + t108, t108, t6, qJDD(2) * mrSges(4,2) - t93 * mrSges(4,3) - t100, t98, t94, -t99;];
f_new  = t3;
