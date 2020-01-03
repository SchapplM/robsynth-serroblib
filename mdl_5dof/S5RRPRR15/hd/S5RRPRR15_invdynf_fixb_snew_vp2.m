% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR15
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
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
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
% Datum: 2019-12-31 20:43
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR15_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR15_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:41:14
% EndTime: 2019-12-31 20:41:17
% DurationCPUTime: 1.08s
% Computational Cost: add. (9569->163), mult. (19674->205), div. (0->0), fcn. (11365->8), ass. (0->82)
t121 = -2 * qJD(3);
t81 = sin(qJ(2));
t73 = t81 * qJD(1);
t66 = mrSges(4,1) * t73 + qJD(2) * mrSges(4,2);
t108 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t73 - t66;
t113 = mrSges(3,1) - mrSges(4,2);
t88 = qJD(1) ^ 2;
t115 = t88 * pkin(6);
t106 = qJD(1) * qJD(2);
t85 = cos(qJ(2));
t101 = t85 * t106;
t60 = t81 * qJDD(1) + t101;
t102 = t81 * t106;
t61 = t85 * qJDD(1) - t102;
t107 = qJD(1) * t85;
t64 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t107;
t82 = sin(qJ(1));
t86 = cos(qJ(1));
t104 = t82 * g(1) - t86 * g(2);
t97 = -qJDD(1) * pkin(1) - t104;
t67 = pkin(3) * t73 - qJD(2) * pkin(7);
t78 = t85 ^ 2;
t91 = pkin(2) * t102 + t73 * t121 + (-t60 - t101) * qJ(3) + t97;
t20 = -t67 * t73 + (-pkin(3) * t78 - pkin(6)) * t88 + (-pkin(2) - pkin(7)) * t61 + t91;
t99 = -t86 * g(1) - t82 * g(2);
t49 = -t88 * pkin(1) + qJDD(1) * pkin(6) + t99;
t110 = -t85 * g(3) - t81 * t49;
t57 = (-pkin(2) * t85 - qJ(3) * t81) * qJD(1);
t87 = qJD(2) ^ 2;
t31 = -qJDD(2) * pkin(2) - t87 * qJ(3) + t57 * t73 + qJDD(3) - t110;
t25 = (-t81 * t85 * t88 - qJDD(2)) * pkin(7) + (t60 - t101) * pkin(3) + t31;
t80 = sin(qJ(4));
t84 = cos(qJ(4));
t100 = -t80 * t20 + t84 * t25;
t55 = -t80 * qJD(2) - t84 * t107;
t37 = t55 * qJD(4) + t84 * qJDD(2) - t80 * t61;
t54 = qJDD(4) + t60;
t56 = t84 * qJD(2) - t80 * t107;
t70 = t73 + qJD(4);
t12 = (t55 * t70 - t37) * pkin(8) + (t55 * t56 + t54) * pkin(4) + t100;
t111 = t84 * t20 + t80 * t25;
t36 = -t56 * qJD(4) - t80 * qJDD(2) - t84 * t61;
t43 = t70 * pkin(4) - t56 * pkin(8);
t53 = t55 ^ 2;
t13 = -t53 * pkin(4) + t36 * pkin(8) - t70 * t43 + t111;
t79 = sin(qJ(5));
t83 = cos(qJ(5));
t38 = t83 * t55 - t79 * t56;
t18 = t38 * qJD(5) + t79 * t36 + t83 * t37;
t39 = t79 * t55 + t83 * t56;
t27 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t68 = qJD(5) + t70;
t32 = -t68 * mrSges(6,2) + t38 * mrSges(6,3);
t50 = qJDD(5) + t54;
t10 = m(6) * (t83 * t12 - t79 * t13) - t18 * mrSges(6,3) + t50 * mrSges(6,1) - t39 * t27 + t68 * t32;
t17 = -t39 * qJD(5) + t83 * t36 - t79 * t37;
t33 = t68 * mrSges(6,1) - t39 * mrSges(6,3);
t11 = m(6) * (t79 * t12 + t83 * t13) + t17 * mrSges(6,3) - t50 * mrSges(6,2) + t38 * t27 - t68 * t33;
t40 = -t55 * mrSges(5,1) + t56 * mrSges(5,2);
t41 = -t70 * mrSges(5,2) + t55 * mrSges(5,3);
t6 = m(5) * t100 + t54 * mrSges(5,1) - t37 * mrSges(5,3) + t83 * t10 + t79 * t11 - t56 * t40 + t70 * t41;
t65 = -mrSges(4,1) * t107 - qJD(2) * mrSges(4,3);
t42 = t70 * mrSges(5,1) - t56 * mrSges(5,3);
t7 = m(5) * t111 - t54 * mrSges(5,2) + t36 * mrSges(5,3) - t79 * t10 + t83 * t11 + t55 * t40 - t70 * t42;
t98 = t80 * t6 - t84 * t7 - m(4) * (-t61 * pkin(2) - t115 + t91) - t65 * t107 + t60 * mrSges(4,3);
t120 = (t108 * t81 - t85 * t64) * qJD(1) - t113 * t61 + m(3) * (t97 - t115) + t60 * mrSges(3,2) - t98;
t103 = -t81 * g(3) + t85 * t49;
t119 = t87 * pkin(2) - qJDD(2) * qJ(3) + qJD(2) * t121 - t57 * t107 - t103;
t58 = (mrSges(4,2) * t85 - mrSges(4,3) * t81) * qJD(1);
t109 = t58 + (-mrSges(3,1) * t85 + mrSges(3,2) * t81) * qJD(1);
t112 = mrSges(3,3) + mrSges(4,1);
t96 = -m(4) * t31 - t84 * t6 - t80 * t7;
t4 = m(3) * t110 - t112 * t60 + t113 * qJDD(2) + (t64 - t65) * qJD(2) - t109 * t73 + t96;
t92 = -t78 * t88 * pkin(7) + t61 * pkin(3) + qJD(2) * t67 - t119;
t94 = -t17 * mrSges(6,1) - t38 * t32 + m(6) * (-t36 * pkin(4) - t53 * pkin(8) + t56 * t43 + t92) + t18 * mrSges(6,2) + t39 * t33;
t90 = m(5) * t92 - t36 * mrSges(5,1) + t37 * mrSges(5,2) - t55 * t41 + t56 * t42 + t94;
t89 = -m(4) * t119 + t90;
t9 = t112 * t61 + (-mrSges(3,2) + mrSges(4,3)) * qJDD(2) - t108 * qJD(2) + m(3) * t103 + t109 * t107 + t89;
t116 = t85 * t4 + t81 * t9;
t2 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t88 * mrSges(2,2) - t120;
t1 = m(2) * t99 - t88 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t81 * t4 + t85 * t9;
t3 = [-m(1) * g(1) + t86 * t1 - t82 * t2, t1, t9, t61 * mrSges(4,2) - t66 * t73 - t98, t7, t11; -m(1) * g(2) + t82 * t1 + t86 * t2, t2, t4, -t61 * mrSges(4,1) - qJDD(2) * mrSges(4,3) - qJD(2) * t66 - t58 * t107 - t89, t6, t10; (-m(1) - m(2)) * g(3) + t116, -m(2) * g(3) + t116, t120, t60 * mrSges(4,1) + qJDD(2) * mrSges(4,2) + qJD(2) * t65 + t58 * t73 - t96, t90, t94;];
f_new = t3;
