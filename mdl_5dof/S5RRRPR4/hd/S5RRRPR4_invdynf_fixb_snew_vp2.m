% Calculate vector of cutting forces with Newton-Euler
% S5RRRPR4
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
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:12
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:11:02
% EndTime: 2019-12-31 21:11:04
% DurationCPUTime: 0.74s
% Computational Cost: add. (7178->141), mult. (9140->181), div. (0->0), fcn. (4687->8), ass. (0->71)
t64 = qJD(1) + qJD(2);
t60 = t64 ^ 2;
t100 = t60 * pkin(7);
t62 = qJDD(1) + qJDD(2);
t69 = sin(qJ(3));
t73 = cos(qJ(3));
t90 = qJD(3) * t73;
t89 = t64 * t90;
t42 = t69 * t62 + t89;
t99 = t64 * t69;
t43 = -qJD(3) * t99 + t73 * t62;
t51 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t99;
t52 = -qJD(3) * mrSges(5,1) + mrSges(5,2) * t99;
t71 = sin(qJ(1));
t75 = cos(qJ(1));
t88 = t71 * g(1) - t75 * g(2);
t49 = qJDD(1) * pkin(1) + t88;
t77 = qJD(1) ^ 2;
t84 = -t75 * g(1) - t71 * g(2);
t50 = -t77 * pkin(1) + t84;
t70 = sin(qJ(2));
t74 = cos(qJ(2));
t94 = t74 * t49 - t70 * t50;
t86 = t62 * pkin(2) + t94;
t82 = -t42 * qJ(4) - t86;
t103 = 2 * qJD(4);
t68 = sin(qJ(5));
t72 = cos(qJ(5));
t37 = (-t68 * t73 + t69 * t72) * t64;
t22 = -t37 * qJD(5) - t68 * t42 - t72 * t43;
t36 = (-t68 * t69 - t72 * t73) * t64;
t23 = t36 * qJD(5) + t72 * t42 - t68 * t43;
t63 = -qJD(3) + qJD(5);
t31 = -t63 * mrSges(6,2) + t36 * mrSges(6,3);
t32 = t63 * mrSges(6,1) - t37 * mrSges(6,3);
t55 = -qJD(3) * pkin(4) - pkin(8) * t99;
t67 = t73 ^ 2;
t85 = m(6) * ((-pkin(8) * t67 + pkin(7)) * t60 + (pkin(3) + pkin(4)) * t43 + (qJ(4) * t90 + (-pkin(3) * qJD(3) + t103 + t55) * t69) * t64 - t82) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t37 * t32 - t36 * t31;
t83 = m(5) * (-t43 * pkin(3) - t100 + (-0.2e1 * qJD(4) * t69 + (pkin(3) * t69 - qJ(4) * t73) * qJD(3)) * t64 + t82) - t43 * mrSges(5,1) - t85;
t98 = t64 * t73;
t54 = mrSges(5,2) * t98 + qJD(3) * mrSges(5,3);
t91 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t98 + t54;
t106 = (-t91 * t73 + (t51 - t52) * t69) * t64 + (mrSges(4,2) - mrSges(5,3)) * t42 + m(4) * (-t86 - t100) - t43 * mrSges(4,1) + t83;
t102 = -m(2) - m(3);
t41 = (-mrSges(4,1) * t73 + mrSges(4,2) * t69) * t64;
t39 = (-pkin(3) * t73 - qJ(4) * t69) * t64;
t76 = qJD(3) ^ 2;
t93 = t70 * t49 + t74 * t50;
t30 = -t60 * pkin(2) + t62 * pkin(7) + t93;
t87 = -t69 * g(3) + t73 * t30;
t79 = -t76 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t103 + t39 * t98 + t87;
t14 = -t67 * t60 * pkin(4) - t43 * pkin(8) + qJD(3) * t55 + t79;
t95 = -t73 * g(3) - t69 * t30;
t19 = -qJDD(3) * pkin(3) - t76 * qJ(4) + t39 * t99 + qJDD(4) - t95;
t15 = (-t42 + t89) * pkin(8) + (-t60 * t69 * t73 - qJDD(3)) * pkin(4) + t19;
t28 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t61 = -qJDD(3) + qJDD(5);
t10 = m(6) * (-t68 * t14 + t72 * t15) - t23 * mrSges(6,3) + t61 * mrSges(6,1) - t37 * t28 + t63 * t31;
t11 = m(6) * (t72 * t14 + t68 * t15) + t22 * mrSges(6,3) - t61 * mrSges(6,2) + t36 * t28 - t63 * t32;
t40 = (-mrSges(5,1) * t73 - mrSges(5,3) * t69) * t64;
t81 = m(5) * t79 + qJDD(3) * mrSges(5,3) + qJD(3) * t52 - t68 * t10 + t72 * t11 + t40 * t98;
t96 = mrSges(4,3) + mrSges(5,2);
t6 = m(4) * t87 - qJDD(3) * mrSges(4,2) - qJD(3) * t51 + t41 * t98 + t96 * t43 + t81;
t80 = -m(5) * t19 - t72 * t10 - t68 * t11;
t7 = m(4) * t95 + (-t40 - t41) * t99 - t96 * t42 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + t91 * qJD(3) + t80;
t101 = t69 * t6 + t73 * t7;
t8 = m(3) * t94 + t62 * mrSges(3,1) - t60 * mrSges(3,2) - t106;
t3 = m(3) * t93 - t60 * mrSges(3,1) - t62 * mrSges(3,2) + t73 * t6 - t69 * t7;
t2 = m(2) * t84 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t74 * t3 - t70 * t8;
t1 = m(2) * t88 + qJDD(1) * mrSges(2,1) - t77 * mrSges(2,2) + t70 * t3 + t74 * t8;
t4 = [-m(1) * g(1) - t71 * t1 + t75 * t2, t2, t3, t6, t43 * mrSges(5,2) + t81, t11; -m(1) * g(2) + t75 * t1 + t71 * t2, t1, t8, t7, -t42 * mrSges(5,3) + (-t69 * t52 - t73 * t54) * t64 + t83, t10; (-m(1) + t102) * g(3) + t101, t102 * g(3) + t101, -m(3) * g(3) + t101, t106, -qJDD(3) * mrSges(5,1) + t42 * mrSges(5,2) - qJD(3) * t54 + t40 * t99 - t80, t85;];
f_new = t4;
