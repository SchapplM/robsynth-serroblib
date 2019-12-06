% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR4
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2019-12-05 17:45
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:44:13
% EndTime: 2019-12-05 17:44:19
% DurationCPUTime: 2.16s
% Computational Cost: add. (17152->147), mult. (46826->210), div. (0->0), fcn. (31990->10), ass. (0->89)
t77 = sin(pkin(8));
t112 = qJD(1) * t77;
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t102 = t85 * g(2) + t82 * g(3);
t86 = qJD(1) ^ 2;
t90 = -t86 * qJ(2) + qJDD(2) - t102;
t79 = cos(pkin(8));
t99 = -pkin(2) * t79 - qJ(3) * t77;
t121 = -0.2e1 * qJD(3) * t112 + (-pkin(1) + t99) * qJDD(1) + t90;
t74 = t77 ^ 2;
t120 = (t79 ^ 2 + t74) * mrSges(3,3);
t76 = sin(pkin(9));
t78 = cos(pkin(9));
t100 = mrSges(4,1) * t76 + mrSges(4,2) * t78;
t106 = t82 * g(2) - t85 * g(3);
t58 = -t86 * pkin(1) + qJDD(1) * qJ(2) + t106;
t114 = -t79 * g(1) - t77 * t58;
t101 = -t79 * mrSges(3,1) + t77 * mrSges(3,2);
t61 = t101 * qJD(1);
t108 = 0.2e1 * qJD(1) * qJD(2);
t60 = t99 * qJD(1);
t35 = t77 * t108 + t60 * t112 + qJDD(3) - t114;
t81 = sin(qJ(4));
t84 = cos(qJ(4));
t93 = (-t76 * t81 + t78 * t84) * t77;
t50 = qJD(1) * t93;
t94 = (-t76 * t84 - t78 * t81) * t77;
t39 = -t50 * qJD(4) + qJDD(1) * t94;
t49 = qJD(1) * t94;
t40 = t49 * qJD(4) + qJDD(1) * t93;
t111 = t79 * qJD(1);
t68 = qJD(4) - t111;
t41 = -t68 * mrSges(5,2) + t49 * mrSges(5,3);
t42 = t68 * mrSges(5,1) - t50 * mrSges(5,3);
t118 = t76 * t77;
t107 = qJDD(1) * t118;
t110 = t76 ^ 2 * t74 * t86;
t117 = t77 * t78;
t97 = -pkin(3) * t79 - pkin(6) * t117;
t59 = t97 * qJD(1);
t89 = t78 * t59 * t112 + pkin(3) * t107 - pkin(6) * t110 + t35;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t34 = t80 * t49 + t83 * t50;
t18 = -t34 * qJD(5) + t83 * t39 - t80 * t40;
t33 = t83 * t49 - t80 * t50;
t19 = t33 * qJD(5) + t80 * t39 + t83 * t40;
t66 = qJD(5) + t68;
t29 = -t66 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = t66 * mrSges(6,1) - t34 * mrSges(6,3);
t45 = t68 * pkin(4) - t50 * pkin(7);
t48 = t49 ^ 2;
t91 = t18 * mrSges(6,1) + t33 * t29 - m(6) * (-t39 * pkin(4) - t48 * pkin(7) + t50 * t45 + t89) - t19 * mrSges(6,2) - t34 * t30;
t88 = -m(5) * t89 + t39 * mrSges(5,1) - t40 * mrSges(5,2) + t49 * t41 - t50 * t42 + t91;
t87 = m(4) * t35 - t88;
t95 = t79 * mrSges(4,2) - mrSges(4,3) * t118;
t56 = t95 * qJD(1);
t96 = -t79 * mrSges(4,1) - mrSges(4,3) * t117;
t57 = t96 * qJD(1);
t98 = t56 * t76 + t57 * t78;
t10 = -t87 + ((-mrSges(3,3) - t100) * qJDD(1) + (-0.2e1 * m(3) * qJD(2) - t61 - t98) * qJD(1)) * t77 + m(3) * t114;
t103 = -t77 * g(1) + (t108 + t58) * t79;
t115 = t121 * t78;
t36 = t60 * t111 + t103;
t53 = t100 * t112;
t23 = t97 * qJDD(1) + (-t36 + (-pkin(3) * t74 * t78 + pkin(6) * t77 * t79) * t86) * t76 + t115;
t109 = t121 * t76 + t78 * t36;
t24 = -pkin(3) * t110 - pkin(6) * t107 + t59 * t111 + t109;
t105 = t84 * t23 - t81 * t24;
t67 = -t79 * qJDD(1) + qJDD(4);
t13 = (t49 * t68 - t40) * pkin(7) + (t49 * t50 + t67) * pkin(4) + t105;
t116 = t81 * t23 + t84 * t24;
t14 = -t48 * pkin(4) + t39 * pkin(7) - t68 * t45 + t116;
t26 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t64 = qJDD(5) + t67;
t11 = m(6) * (t83 * t13 - t80 * t14) - t19 * mrSges(6,3) + t64 * mrSges(6,1) - t34 * t26 + t66 * t29;
t12 = m(6) * (t80 * t13 + t83 * t14) + t18 * mrSges(6,3) - t64 * mrSges(6,2) + t33 * t26 - t66 * t30;
t37 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t7 = m(5) * t105 + t67 * mrSges(5,1) - t40 * mrSges(5,3) + t83 * t11 + t80 * t12 - t50 * t37 + t68 * t41;
t8 = m(5) * t116 - t67 * mrSges(5,2) + t39 * mrSges(5,3) - t80 * t11 + t83 * t12 + t49 * t37 - t68 * t42;
t5 = m(4) * (-t76 * t36 + t115) + t81 * t8 + t84 * t7 + t96 * qJDD(1) + (-t53 * t117 - t79 * t56) * qJD(1);
t6 = m(4) * t109 + t84 * t8 - t81 * t7 + t95 * qJDD(1) + (-t53 * t118 + t79 * t57) * qJD(1);
t4 = m(3) * t103 + t78 * t6 - t76 * t5 + (qJDD(1) * mrSges(3,3) + qJD(1) * t61) * t79;
t119 = t79 * t10 + t77 * t4;
t92 = m(3) * (-qJDD(1) * pkin(1) + t90) + t78 * t5 + t76 * t6;
t2 = m(2) * t102 + (-mrSges(2,2) + t120) * t86 + (mrSges(2,1) - t101) * qJDD(1) - t92;
t1 = m(2) * t106 - t86 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t77 * t10 + t79 * t4;
t3 = [(-m(1) - m(2)) * g(1) + t119, t1, t4, t6, t8, t12; -m(1) * g(2) - t82 * t1 - t85 * t2, t2, t10, t5, t7, t11; -m(1) * g(3) + t85 * t1 - t82 * t2, -m(2) * g(1) + t119, t101 * qJDD(1) - t120 * t86 + t92, (t98 * qJD(1) + t100 * qJDD(1)) * t77 + t87, -t88, -t91;];
f_new = t3;
