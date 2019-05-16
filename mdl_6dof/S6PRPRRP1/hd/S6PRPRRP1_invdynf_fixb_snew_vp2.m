% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRP1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-05-04 23:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 23:29:10
% EndTime: 2019-05-04 23:29:14
% DurationCPUTime: 1.49s
% Computational Cost: add. (19823->146), mult. (36337->190), div. (0->0), fcn. (24713->12), ass. (0->81)
t83 = sin(qJ(4));
t108 = qJD(2) * t83;
t82 = sin(qJ(5));
t85 = cos(qJ(5));
t60 = t85 * qJD(4) - t108 * t82;
t106 = qJD(2) * qJD(4);
t86 = cos(qJ(4));
t98 = t86 * t106;
t64 = t83 * qJDD(2) + t98;
t41 = t60 * qJD(5) + t82 * qJDD(4) + t85 * t64;
t107 = t86 * qJD(2);
t72 = qJD(5) - t107;
t47 = -t72 * mrSges(7,2) + t60 * mrSges(7,3);
t99 = t83 * t106;
t65 = t86 * qJDD(2) - t99;
t58 = qJDD(5) - t65;
t61 = t82 * qJD(4) + t108 * t85;
t78 = sin(pkin(6));
t87 = cos(qJ(2));
t114 = t78 * t87;
t77 = sin(pkin(10));
t80 = cos(pkin(10));
t66 = t77 * g(1) - t80 * g(2);
t81 = cos(pkin(6));
t116 = t66 * t81;
t67 = -t80 * g(1) - t77 * g(2);
t75 = -g(3) + qJDD(1);
t84 = sin(qJ(2));
t93 = t75 * t114 + t87 * t116 - t84 * t67;
t35 = qJDD(2) * pkin(2) + t93;
t115 = t78 * t84;
t101 = t75 * t115 + t84 * t116 + t87 * t67;
t89 = qJD(2) ^ 2;
t36 = -t89 * pkin(2) + t101;
t76 = sin(pkin(11));
t79 = cos(pkin(11));
t110 = t76 * t35 + t79 * t36;
t30 = -t89 * pkin(3) + qJDD(2) * pkin(8) + t110;
t94 = -t78 * t66 + t81 * t75;
t52 = qJDD(3) + t94;
t111 = t86 * t30 + t83 * t52;
t63 = (-pkin(4) * t86 - pkin(9) * t83) * qJD(2);
t88 = qJD(4) ^ 2;
t23 = -t88 * pkin(4) + qJDD(4) * pkin(9) + t107 * t63 + t111;
t95 = t79 * t35 - t76 * t36;
t29 = -qJDD(2) * pkin(3) - t89 * pkin(8) - t95;
t26 = (-t64 - t98) * pkin(9) + (-t65 + t99) * pkin(4) + t29;
t97 = -t82 * t23 + t85 * t26;
t105 = m(7) * (-0.2e1 * qJD(6) * t61 + (t60 * t72 - t41) * qJ(6) + (t60 * t61 + t58) * pkin(5) + t97) + t72 * t47 + t58 * mrSges(7,1);
t43 = -t60 * mrSges(7,1) + t61 * mrSges(7,2);
t44 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t48 = -t72 * mrSges(6,2) + t60 * mrSges(6,3);
t13 = m(6) * t97 + t58 * mrSges(6,1) + t72 * t48 + (-t44 - t43) * t61 + (-mrSges(6,3) - mrSges(7,3)) * t41 + t105;
t112 = t85 * t23 + t82 * t26;
t40 = -t61 * qJD(5) + t85 * qJDD(4) - t82 * t64;
t49 = t72 * pkin(5) - t61 * qJ(6);
t57 = t60 ^ 2;
t104 = m(7) * (-t57 * pkin(5) + t40 * qJ(6) + 0.2e1 * qJD(6) * t60 - t72 * t49 + t112) + t60 * t43 + t40 * mrSges(7,3);
t50 = t72 * mrSges(7,1) - t61 * mrSges(7,3);
t51 = t72 * mrSges(6,1) - t61 * mrSges(6,3);
t14 = m(6) * t112 + t40 * mrSges(6,3) + t60 * t44 + (-t51 - t50) * t72 + (-mrSges(6,2) - mrSges(7,2)) * t58 + t104;
t68 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108;
t69 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t107;
t118 = m(5) * t29 - t65 * mrSges(5,1) + t64 * mrSges(5,2) + t85 * t13 + t82 * t14 + (t83 * t68 - t86 * t69) * qJD(2);
t96 = -t83 * t30 + t86 * t52;
t22 = -qJDD(4) * pkin(4) - t88 * pkin(9) + t63 * t108 - t96;
t103 = m(7) * (-t40 * pkin(5) - t57 * qJ(6) + t61 * t49 + qJDD(6) + t22) + t41 * mrSges(7,2) + t61 * t50;
t117 = m(6) * t22 + t41 * mrSges(6,2) - (t48 + t47) * t60 - (mrSges(6,1) + mrSges(7,1)) * t40 + t61 * t51 + t103;
t62 = (-mrSges(5,1) * t86 + mrSges(5,2) * t83) * qJD(2);
t12 = m(5) * t111 - qJDD(4) * mrSges(5,2) + t65 * mrSges(5,3) - qJD(4) * t68 + t107 * t62 - t82 * t13 + t85 * t14;
t16 = m(5) * t96 + qJDD(4) * mrSges(5,1) - t64 * mrSges(5,3) + qJD(4) * t69 - t108 * t62 - t117;
t102 = m(4) * t52 + t83 * t12 + t86 * t16;
t10 = m(4) * t95 + qJDD(2) * mrSges(4,1) - t89 * mrSges(4,2) - t118;
t7 = m(4) * t110 - t89 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t86 * t12 - t83 * t16;
t5 = m(3) * t93 + qJDD(2) * mrSges(3,1) - t89 * mrSges(3,2) + t79 * t10 + t76 * t7;
t6 = m(3) * t101 - t89 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t76 * t10 + t79 * t7;
t9 = m(3) * t94 + t102;
t100 = m(2) * t75 + t5 * t114 + t6 * t115 + t81 * t9;
t2 = m(2) * t67 - t84 * t5 + t87 * t6;
t1 = m(2) * t66 - t78 * t9 + (t5 * t87 + t6 * t84) * t81;
t3 = [-m(1) * g(1) - t77 * t1 + t80 * t2, t2, t6, t7, t12, t14, -t58 * mrSges(7,2) - t72 * t50 + t104; -m(1) * g(2) + t80 * t1 + t77 * t2, t1, t5, t10, t16, t13, -t41 * mrSges(7,3) - t61 * t43 + t105; -m(1) * g(3) + t100, t100, t9, t102, t118, t117, -t40 * mrSges(7,1) - t60 * t47 + t103;];
f_new  = t3;
