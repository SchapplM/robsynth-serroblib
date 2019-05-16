% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:09
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:05:08
% EndTime: 2019-05-05 00:05:13
% DurationCPUTime: 2.62s
% Computational Cost: add. (38690->151), mult. (72811->209), div. (0->0), fcn. (51996->14), ass. (0->90)
t89 = sin(qJ(5));
t90 = sin(qJ(4));
t93 = cos(qJ(5));
t94 = cos(qJ(4));
t61 = (t89 * t90 - t93 * t94) * qJD(2);
t84 = sin(pkin(6));
t95 = cos(qJ(2));
t119 = t84 * t95;
t83 = sin(pkin(11));
t86 = cos(pkin(11));
t68 = t83 * g(1) - t86 * g(2);
t87 = cos(pkin(6));
t121 = t68 * t87;
t69 = -t86 * g(1) - t83 * g(2);
t81 = -g(3) + qJDD(1);
t91 = sin(qJ(2));
t105 = t81 * t119 + t95 * t121 - t91 * t69;
t45 = qJDD(2) * pkin(2) + t105;
t120 = t84 * t91;
t111 = t81 * t120 + t91 * t121 + t95 * t69;
t96 = qJD(2) ^ 2;
t46 = -t96 * pkin(2) + t111;
t82 = sin(pkin(12));
t85 = cos(pkin(12));
t116 = t82 * t45 + t85 * t46;
t32 = -t96 * pkin(3) + qJDD(2) * pkin(8) + t116;
t106 = -t84 * t68 + t87 * t81;
t57 = qJDD(3) + t106;
t108 = -t90 * t32 + t94 * t57;
t113 = qJD(2) * qJD(4);
t109 = t94 * t113;
t66 = t90 * qJDD(2) + t109;
t26 = (-t66 + t109) * pkin(9) + (t90 * t94 * t96 + qJDD(4)) * pkin(4) + t108;
t117 = t94 * t32 + t90 * t57;
t67 = t94 * qJDD(2) - t90 * t113;
t115 = qJD(2) * t90;
t74 = qJD(4) * pkin(4) - pkin(9) * t115;
t80 = t94 ^ 2;
t27 = -t80 * t96 * pkin(4) + t67 * pkin(9) - qJD(4) * t74 + t117;
t118 = t89 * t26 + t93 * t27;
t62 = (t89 * t94 + t90 * t93) * qJD(2);
t49 = t61 * pkin(5) - t62 * pkin(10);
t79 = qJD(4) + qJD(5);
t77 = t79 ^ 2;
t78 = qJDD(4) + qJDD(5);
t22 = -t77 * pkin(5) + t78 * pkin(10) - t61 * t49 + t118;
t38 = -t62 * qJD(5) - t89 * t66 + t93 * t67;
t39 = -t61 * qJD(5) + t93 * t66 + t89 * t67;
t107 = t85 * t45 - t82 * t46;
t101 = -qJDD(2) * pkin(3) - t107;
t98 = -t67 * pkin(4) + t74 * t115 + (-pkin(9) * t80 - pkin(8)) * t96 + t101;
t23 = (t61 * t79 - t39) * pkin(10) + (t62 * t79 - t38) * pkin(5) + t98;
t88 = sin(qJ(6));
t92 = cos(qJ(6));
t50 = -t88 * t62 + t92 * t79;
t34 = t50 * qJD(6) + t92 * t39 + t88 * t78;
t51 = t92 * t62 + t88 * t79;
t35 = -t50 * mrSges(7,1) + t51 * mrSges(7,2);
t37 = qJDD(6) - t38;
t58 = qJD(6) + t61;
t40 = -t58 * mrSges(7,2) + t50 * mrSges(7,3);
t19 = m(7) * (-t88 * t22 + t92 * t23) - t34 * mrSges(7,3) + t37 * mrSges(7,1) - t51 * t35 + t58 * t40;
t33 = -t51 * qJD(6) - t88 * t39 + t92 * t78;
t41 = t58 * mrSges(7,1) - t51 * mrSges(7,3);
t20 = m(7) * (t92 * t22 + t88 * t23) + t33 * mrSges(7,3) - t37 * mrSges(7,2) + t50 * t35 - t58 * t41;
t55 = -t79 * mrSges(6,2) - t61 * mrSges(6,3);
t56 = t79 * mrSges(6,1) - t62 * mrSges(6,3);
t100 = -m(6) * t98 + t38 * mrSges(6,1) - t39 * mrSges(6,2) - t92 * t19 - t88 * t20 - t61 * t55 - t62 * t56;
t70 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t115;
t114 = qJD(2) * t94;
t71 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t114;
t122 = (t90 * t70 - t94 * t71) * qJD(2) + m(5) * (-t96 * pkin(8) + t101) - t67 * mrSges(5,1) + t66 * mrSges(5,2) - t100;
t48 = t61 * mrSges(6,1) + t62 * mrSges(6,2);
t15 = m(6) * t118 - t78 * mrSges(6,2) + t38 * mrSges(6,3) - t88 * t19 + t92 * t20 - t61 * t48 - t79 * t56;
t104 = t93 * t26 - t89 * t27;
t99 = m(7) * (-t78 * pkin(5) - t77 * pkin(10) + t62 * t49 - t104) - t33 * mrSges(7,1) + t34 * mrSges(7,2) - t50 * t40 + t51 * t41;
t16 = m(6) * t104 + t78 * mrSges(6,1) - t39 * mrSges(6,3) - t62 * t48 + t79 * t55 - t99;
t65 = (-mrSges(5,1) * t94 + mrSges(5,2) * t90) * qJD(2);
t12 = m(5) * t108 + qJDD(4) * mrSges(5,1) - t66 * mrSges(5,3) + qJD(4) * t71 - t65 * t115 + t89 * t15 + t93 * t16;
t13 = m(5) * t117 - qJDD(4) * mrSges(5,2) + t67 * mrSges(5,3) - qJD(4) * t70 + t65 * t114 + t93 * t15 - t89 * t16;
t112 = m(4) * t57 + t94 * t12 + t90 * t13;
t14 = m(4) * t107 + qJDD(2) * mrSges(4,1) - t96 * mrSges(4,2) - t122;
t7 = m(4) * t116 - t96 * mrSges(4,1) - qJDD(2) * mrSges(4,2) - t90 * t12 + t94 * t13;
t5 = m(3) * t105 + qJDD(2) * mrSges(3,1) - t96 * mrSges(3,2) + t85 * t14 + t82 * t7;
t6 = m(3) * t111 - t96 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t82 * t14 + t85 * t7;
t9 = m(3) * t106 + t112;
t110 = m(2) * t81 + t5 * t119 + t6 * t120 + t87 * t9;
t2 = m(2) * t69 - t91 * t5 + t95 * t6;
t1 = m(2) * t68 - t84 * t9 + (t5 * t95 + t6 * t91) * t87;
t3 = [-m(1) * g(1) - t83 * t1 + t86 * t2, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t86 * t1 + t83 * t2, t1, t5, t14, t12, t16, t19; -m(1) * g(3) + t110, t110, t9, t112, t122, -t100, t99;];
f_new  = t3;
