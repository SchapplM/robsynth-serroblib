% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP5
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1]';
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
% Datum: 2019-05-04 23:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP5_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:54:15
% EndTime: 2019-05-04 23:54:17
% DurationCPUTime: 0.91s
% Computational Cost: add. (10384->146), mult. (19106->182), div. (0->0), fcn. (11926->10), ass. (0->79)
t77 = sin(pkin(6));
t82 = sin(qJ(2));
t116 = t77 * t82;
t76 = sin(pkin(10));
t78 = cos(pkin(10));
t66 = t76 * g(1) - t78 * g(2);
t79 = cos(pkin(6));
t117 = t66 * t79;
t67 = -t78 * g(1) - t76 * g(2);
t75 = -g(3) + qJDD(1);
t85 = cos(qJ(2));
t102 = t75 * t116 + t82 * t117 + t85 * t67;
t122 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t102;
t121 = -t82 * t67 + (t75 * t77 + t117) * t85;
t84 = cos(qJ(4));
t109 = qJD(2) * t84;
t81 = sin(qJ(4));
t63 = (pkin(4) * t81 - pkin(9) * t84) * qJD(2);
t86 = qJD(4) ^ 2;
t119 = -pkin(2) - pkin(8);
t87 = qJD(2) ^ 2;
t88 = -t87 * qJ(3) + qJDD(3) - t121;
t29 = t119 * qJDD(2) + t88;
t49 = -t77 * t66 + t79 * t75;
t97 = t84 * t29 - t81 * t49;
t21 = -qJDD(4) * pkin(4) - t86 * pkin(9) + t63 * t109 - t97;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t61 = t80 * qJD(4) + t83 * t109;
t107 = qJD(2) * qJD(4);
t100 = t81 * t107;
t65 = t84 * qJDD(2) - t100;
t36 = -t61 * qJD(5) + t83 * qJDD(4) - t80 * t65;
t60 = t83 * qJD(4) - t80 * t109;
t37 = t60 * qJD(5) + t80 * qJDD(4) + t83 * t65;
t108 = t81 * qJD(2);
t71 = qJD(5) + t108;
t45 = t71 * pkin(5) - t61 * qJ(6);
t46 = t71 * mrSges(7,1) - t61 * mrSges(7,3);
t56 = t60 ^ 2;
t103 = m(7) * (-t36 * pkin(5) - t56 * qJ(6) + t61 * t45 + qJDD(6) + t21) + t37 * mrSges(7,2) + t61 * t46;
t43 = -t71 * mrSges(7,2) + t60 * mrSges(7,3);
t44 = -t71 * mrSges(6,2) + t60 * mrSges(6,3);
t47 = t71 * mrSges(6,1) - t61 * mrSges(6,3);
t120 = m(6) * t21 + t37 * mrSges(6,2) - (t44 + t43) * t60 - (mrSges(6,1) + mrSges(7,1)) * t36 + t61 * t47 + t103;
t113 = (-mrSges(3,2) + mrSges(4,3));
t115 = mrSges(3,1) - mrSges(4,2);
t111 = t81 * t29 + t84 * t49;
t99 = t84 * t107;
t64 = -t81 * qJDD(2) - t99;
t57 = qJDD(5) - t64;
t22 = -t86 * pkin(4) + qJDD(4) * pkin(9) - t63 * t108 + t111;
t91 = t119 * t87 - t122;
t25 = (-t65 + t100) * pkin(9) + (-t64 + t99) * pkin(4) + t91;
t98 = -t80 * t22 + t83 * t25;
t105 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t71 - t37) * qJ(6) + (t60 * t61 + t57) * pkin(5) + t98) + t71 * t43 + t57 * mrSges(7,1);
t39 = -t60 * mrSges(7,1) + t61 * mrSges(7,2);
t40 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t12 = m(6) * t98 + t57 * mrSges(6,1) + t71 * t44 + (-t40 - t39) * t61 + (-mrSges(6,3) - mrSges(7,3)) * t37 + t105;
t112 = t83 * t22 + t80 * t25;
t104 = m(7) * (-t56 * pkin(5) + t36 * qJ(6) + 0.2e1 * qJD(6) * t60 - t71 * t45 + t112) + t60 * t39 + t36 * mrSges(7,3);
t14 = m(6) * t112 + t36 * mrSges(6,3) + t60 * t40 + (-t47 - t46) * t71 + (-mrSges(6,2) - mrSges(7,2)) * t57 + t104;
t62 = (mrSges(5,1) * t81 + mrSges(5,2) * t84) * qJD(2);
t69 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t109;
t10 = m(5) * t111 - qJDD(4) * mrSges(5,2) + t64 * mrSges(5,3) - qJD(4) * t69 - t62 * t108 - t80 * t12 + t83 * t14;
t68 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t108;
t15 = m(5) * t97 + qJDD(4) * mrSges(5,1) - t65 * mrSges(5,3) + qJD(4) * t68 - t62 * t109 - t120;
t93 = -m(4) * (-qJDD(2) * pkin(2) + t88) - t81 * t10 - t84 * t15;
t4 = m(3) * t121 + t115 * qJDD(2) + (t113 * t87) + t93;
t118 = t4 * t85;
t95 = m(4) * t49 + t84 * t10 - t81 * t15;
t6 = m(3) * t49 + t95;
t92 = m(5) * t91 - t64 * mrSges(5,1) + t65 * mrSges(5,2) + t68 * t108 + t69 * t109 + t83 * t12 + t80 * t14;
t89 = -m(4) * (t87 * pkin(2) + t122) + t92;
t8 = m(3) * t102 + t113 * qJDD(2) - t115 * t87 + t89;
t101 = m(2) * t75 + t8 * t116 + t77 * t118 + t79 * t6;
t2 = m(2) * t67 - t82 * t4 + t85 * t8;
t1 = m(2) * t66 - t77 * t6 + (t8 * t82 + t118) * t79;
t3 = [-m(1) * g(1) - t76 * t1 + t78 * t2, t2, t8, t95, t10, t14, -t57 * mrSges(7,2) - t71 * t46 + t104; -m(1) * g(2) + t78 * t1 + t76 * t2, t1, t4, -(t87 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t89, t15, t12, -t37 * mrSges(7,3) - t61 * t39 + t105; -m(1) * g(3) + t101, t101, t6, qJDD(2) * mrSges(4,2) - t87 * mrSges(4,3) - t93, t92, t120, -t36 * mrSges(7,1) - t60 * t43 + t103;];
f_new  = t3;
