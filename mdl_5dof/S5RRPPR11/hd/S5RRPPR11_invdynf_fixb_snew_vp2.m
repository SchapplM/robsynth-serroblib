% Calculate vector of cutting forces with Newton-Euler
% S5RRPPR11
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
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta4]';
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
% Datum: 2019-12-31 19:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:46:29
% EndTime: 2019-12-31 19:46:31
% DurationCPUTime: 1.03s
% Computational Cost: add. (8749->162), mult. (19132->207), div. (0->0), fcn. (10858->8), ass. (0->80)
t121 = -2 * qJD(3);
t81 = sin(qJ(2));
t107 = t81 * qJD(1);
t68 = mrSges(4,1) * t107 + qJD(2) * mrSges(4,2);
t109 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t107 - t68;
t113 = mrSges(3,1) - mrSges(4,2);
t87 = qJD(1) ^ 2;
t115 = t87 * pkin(6);
t106 = qJD(1) * qJD(2);
t84 = cos(qJ(2));
t100 = t84 * t106;
t61 = t81 * qJDD(1) + t100;
t101 = t81 * t106;
t62 = t84 * qJDD(1) - t101;
t108 = qJD(1) * t84;
t65 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t108;
t82 = sin(qJ(1));
t85 = cos(qJ(1));
t103 = t82 * g(1) - t85 * g(2);
t96 = -qJDD(1) * pkin(1) - t103;
t78 = sin(pkin(8));
t79 = cos(pkin(8));
t42 = t79 * qJDD(2) - t78 * t62;
t54 = -t78 * qJD(2) - t79 * t108;
t55 = t79 * qJD(2) - t78 * t108;
t66 = pkin(3) * t107 - qJD(2) * qJ(4);
t77 = t84 ^ 2;
t91 = pkin(2) * t101 + t107 * t121 + (-t61 - t100) * qJ(3) + t96;
t17 = -t66 * t107 + (-pkin(3) * t77 - pkin(6)) * t87 + (-pkin(2) - qJ(4)) * t62 + t91;
t98 = -t85 * g(1) - t82 * g(2);
t51 = -t87 * pkin(1) + qJDD(1) * pkin(6) + t98;
t111 = -t84 * g(3) - t81 * t51;
t58 = (-pkin(2) * t84 - qJ(3) * t81) * qJD(1);
t86 = qJD(2) ^ 2;
t31 = -qJDD(2) * pkin(2) - t86 * qJ(3) + t58 * t107 + qJDD(3) - t111;
t25 = (-t81 * t84 * t87 - qJDD(2)) * qJ(4) + (t61 - t100) * pkin(3) + t31;
t99 = -0.2e1 * qJD(4) * t55 - t78 * t17 + t79 * t25;
t12 = (t54 * t107 - t42) * pkin(7) + (t54 * t55 + t61) * pkin(4) + t99;
t104 = 0.2e1 * qJD(4) * t54 + t79 * t17 + t78 * t25;
t41 = -t78 * qJDD(2) - t79 * t62;
t43 = pkin(4) * t107 - t55 * pkin(7);
t53 = t54 ^ 2;
t13 = -t53 * pkin(4) + t41 * pkin(7) - t43 * t107 + t104;
t80 = sin(qJ(5));
t83 = cos(qJ(5));
t35 = t83 * t54 - t80 * t55;
t20 = t35 * qJD(5) + t80 * t41 + t83 * t42;
t36 = t80 * t54 + t83 * t55;
t27 = -t35 * mrSges(6,1) + t36 * mrSges(6,2);
t70 = qJD(5) + t107;
t32 = -t70 * mrSges(6,2) + t35 * mrSges(6,3);
t57 = qJDD(5) + t61;
t10 = m(6) * (t83 * t12 - t80 * t13) - t20 * mrSges(6,3) + t57 * mrSges(6,1) - t36 * t27 + t70 * t32;
t19 = -t36 * qJD(5) + t83 * t41 - t80 * t42;
t33 = t70 * mrSges(6,1) - t36 * mrSges(6,3);
t11 = m(6) * (t80 * t12 + t83 * t13) + t19 * mrSges(6,3) - t57 * mrSges(6,2) + t35 * t27 - t70 * t33;
t37 = -t54 * mrSges(5,1) + t55 * mrSges(5,2);
t39 = -mrSges(5,2) * t107 + t54 * mrSges(5,3);
t6 = m(5) * t99 + t61 * mrSges(5,1) - t42 * mrSges(5,3) + t83 * t10 + t39 * t107 + t80 * t11 - t55 * t37;
t67 = -mrSges(4,1) * t108 - qJD(2) * mrSges(4,3);
t40 = mrSges(5,1) * t107 - t55 * mrSges(5,3);
t7 = m(5) * t104 - t61 * mrSges(5,2) + t41 * mrSges(5,3) - t80 * t10 - t40 * t107 + t83 * t11 + t54 * t37;
t97 = t78 * t6 - t79 * t7 - m(4) * (-t62 * pkin(2) - t115 + t91) - t67 * t108 + t61 * mrSges(4,3);
t120 = (t109 * t81 - t84 * t65) * qJD(1) - t113 * t62 + m(3) * (t96 - t115) + t61 * mrSges(3,2) - t97;
t102 = -t81 * g(3) + t84 * t51;
t119 = t86 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t121 - t58 * t108 - t102;
t59 = (mrSges(4,2) * t84 - mrSges(4,3) * t81) * qJD(1);
t110 = t59 + (-mrSges(3,1) * t84 + mrSges(3,2) * t81) * qJD(1);
t112 = mrSges(3,3) + mrSges(4,1);
t95 = -m(4) * t31 - t79 * t6 - t78 * t7;
t4 = m(3) * t111 - t112 * t61 + t113 * qJDD(2) + (t65 - t67) * qJD(2) - t110 * t107 + t95;
t89 = -t77 * t87 * qJ(4) + t62 * pkin(3) + qJD(2) * t66 + qJDD(4) - t119;
t93 = -t19 * mrSges(6,1) - t35 * t32 + m(6) * (-t41 * pkin(4) - t53 * pkin(7) + t55 * t43 + t89) + t20 * mrSges(6,2) + t36 * t33;
t90 = m(5) * t89 - t41 * mrSges(5,1) + t42 * mrSges(5,2) - t54 * t39 + t55 * t40 + t93;
t88 = -m(4) * t119 + t90;
t9 = t88 + t112 * t62 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t109 * qJD(2) + m(3) * t102 + t110 * t108;
t116 = t84 * t4 + t81 * t9;
t2 = m(2) * t103 + qJDD(1) * mrSges(2,1) - t87 * mrSges(2,2) - t120;
t1 = m(2) * t98 - t87 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t81 * t4 + t84 * t9;
t3 = [-m(1) * g(1) + t85 * t1 - t82 * t2, t1, t9, t62 * mrSges(4,2) - t68 * t107 - t97, t7, t11; -m(1) * g(2) + t82 * t1 + t85 * t2, t2, t4, -t62 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t68 - t59 * t108 - t88, t6, t10; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t120, t61 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t67 + t59 * t107 - t95, t90, t93;];
f_new = t3;
