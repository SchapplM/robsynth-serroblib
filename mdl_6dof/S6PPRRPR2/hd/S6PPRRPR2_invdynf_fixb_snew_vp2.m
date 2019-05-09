% Calculate vector of cutting forces with Newton-Euler
% S6PPRRPR2
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:19
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:15:45
% EndTime: 2019-05-04 20:15:49
% DurationCPUTime: 1.90s
% Computational Cost: add. (24880->148), mult. (45329->193), div. (0->0), fcn. (32888->14), ass. (0->89)
t71 = sin(pkin(11));
t75 = cos(pkin(11));
t57 = -g(1) * t75 - g(2) * t71;
t70 = sin(pkin(12));
t74 = cos(pkin(12));
t56 = g(1) * t71 - g(2) * t75;
t69 = -g(3) + qJDD(1);
t73 = sin(pkin(6));
t77 = cos(pkin(6));
t95 = t56 * t77 + t69 * t73;
t34 = -t70 * t57 + t95 * t74;
t44 = -t56 * t73 + t69 * t77 + qJDD(2);
t72 = sin(pkin(7));
t76 = cos(pkin(7));
t125 = t34 * t76 + t44 * t72;
t79 = sin(qJ(4));
t104 = qJD(3) * t79;
t82 = cos(qJ(4));
t52 = (mrSges(6,2) * t82 - mrSges(6,3) * t79) * qJD(3);
t106 = t52 + (-mrSges(5,1) * t82 + mrSges(5,2) * t79) * qJD(3);
t108 = mrSges(5,3) + mrSges(6,1);
t109 = mrSges(5,1) - mrSges(6,2);
t30 = -t34 * t72 + t44 * t76;
t113 = t30 * t82;
t35 = t74 * t57 + t95 * t70;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t101 = t125 * t80 + t83 * t35;
t85 = qJD(3) ^ 2;
t28 = -pkin(3) * t85 + qJDD(3) * pkin(9) + t101;
t25 = t79 * t28;
t102 = qJD(3) * qJD(4);
t99 = t82 * t102;
t54 = qJDD(3) * t79 + t99;
t103 = qJD(3) * t82;
t59 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t103;
t60 = -mrSges(6,1) * t103 - qJD(4) * mrSges(6,3);
t116 = pkin(10) * t85;
t118 = -pkin(4) - pkin(10);
t51 = (-pkin(4) * t82 - qJ(5) * t79) * qJD(3);
t84 = qJD(4) ^ 2;
t94 = -qJ(5) * t84 + t51 * t104 + qJDD(5) + t25;
t19 = t54 * pkin(5) + t118 * qJDD(4) + (-pkin(5) * t102 - t79 * t116 - t30) * t82 + t94;
t98 = t79 * t102;
t55 = qJDD(3) * t82 - t98;
t62 = pkin(5) * t104 - qJD(4) * pkin(10);
t68 = t82 ^ 2;
t119 = -2 * qJD(5);
t122 = t125 * t83 - t80 * t35;
t89 = -qJDD(3) * pkin(3) - t122;
t86 = pkin(4) * t98 + t104 * t119 + (-t54 - t99) * qJ(5) + t89;
t22 = -t62 * t104 + (-pkin(5) * t68 - pkin(9)) * t85 + t118 * t55 + t86;
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t49 = -qJD(4) * t78 - t81 * t103;
t39 = t49 * qJD(6) + qJDD(4) * t81 - t55 * t78;
t50 = qJD(4) * t81 - t78 * t103;
t40 = -mrSges(7,1) * t49 + mrSges(7,2) * t50;
t64 = qJD(6) + t104;
t42 = -mrSges(7,2) * t64 + mrSges(7,3) * t49;
t48 = qJDD(6) + t54;
t15 = m(7) * (t19 * t81 - t22 * t78) - t39 * mrSges(7,3) + t48 * mrSges(7,1) - t50 * t40 + t64 * t42;
t38 = -t50 * qJD(6) - qJDD(4) * t78 - t55 * t81;
t43 = mrSges(7,1) * t64 - mrSges(7,3) * t50;
t16 = m(7) * (t19 * t78 + t22 * t81) + t38 * mrSges(7,3) - t48 * mrSges(7,2) + t49 * t40 - t64 * t43;
t92 = -m(6) * (-qJDD(4) * pkin(4) - t113 + t94) - t81 * t15 - t78 * t16;
t12 = m(5) * (-t25 + t113) - t108 * t54 + t109 * qJDD(4) + (t59 - t60) * qJD(4) - t106 * t104 + t92;
t61 = mrSges(6,1) * t104 + qJD(4) * mrSges(6,2);
t105 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t104 - t61;
t107 = t82 * t28 + t79 * t30;
t88 = -pkin(4) * t84 + qJDD(4) * qJ(5) + t51 * t103 + t107;
t91 = -t38 * mrSges(7,1) - t49 * t42 + m(7) * (-t68 * t116 + t55 * pkin(5) + ((2 * qJD(5)) + t62) * qJD(4) + t88) + t39 * mrSges(7,2) + t50 * t43;
t90 = -m(6) * (qJD(4) * t119 - t88) + t91;
t13 = m(5) * t107 + t108 * t55 + (-mrSges(5,2) + mrSges(6,3)) * qJDD(4) - t105 * qJD(4) + t106 * t103 + t90;
t10 = m(4) * t30 + t12 * t82 + t13 * t79;
t117 = pkin(9) * t85;
t93 = t78 * t15 - t81 * t16 - m(6) * (-t55 * pkin(4) - t117 + t86) - t60 * t103 + t54 * mrSges(6,3);
t123 = (t105 * t79 - t82 * t59) * qJD(3) - t109 * t55 + m(5) * (t89 - t117) + t54 * mrSges(5,2) - t93;
t11 = m(4) * t122 + qJDD(3) * mrSges(4,1) - t85 * mrSges(4,2) - t123;
t9 = m(4) * t101 - t85 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t79 * t12 + t82 * t13;
t97 = t11 * t83 + t80 * t9;
t4 = m(3) * t34 - t72 * t10 + t97 * t76;
t8 = m(3) * t35 - t11 * t80 + t83 * t9;
t124 = t4 * t74 + t70 * t8;
t6 = m(3) * t44 + t76 * t10 + t97 * t72;
t100 = m(2) * t69 + t124 * t73 + t77 * t6;
t2 = m(2) * t57 - t4 * t70 + t74 * t8;
t1 = m(2) * t56 + t124 * t77 - t73 * t6;
t3 = [-m(1) * g(1) - t1 * t71 + t2 * t75, t2, t8, t9, t13, t55 * mrSges(6,2) - t61 * t104 - t93, t16; -m(1) * g(2) + t1 * t75 + t2 * t71, t1, t4, t11, t12, -t55 * mrSges(6,1) - qJDD(4) * mrSges(6,3) - qJD(4) * t61 - t52 * t103 - t90, t15; -m(1) * g(3) + t100, t100, t6, t10, t123, t54 * mrSges(6,1) + qJDD(4) * mrSges(6,2) + qJD(4) * t60 + t52 * t104 - t92, t91;];
f_new  = t3;
