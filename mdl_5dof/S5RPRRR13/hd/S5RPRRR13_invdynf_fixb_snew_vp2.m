% Calculate vector of cutting forces with Newton-Euler
% S5RPRRR13
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
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5]';
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
% Datum: 2019-12-31 19:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRRR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:14:35
% EndTime: 2019-12-31 19:14:37
% DurationCPUTime: 0.82s
% Computational Cost: add. (8881->138), mult. (17152->174), div. (0->0), fcn. (10373->8), ass. (0->75)
t69 = sin(qJ(1));
t73 = cos(qJ(1));
t85 = -t73 * g(1) - t69 * g(2);
t99 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t85;
t98 = -m(2) - m(3);
t97 = (-pkin(1) - pkin(6));
t96 = (mrSges(2,1) - mrSges(3,2));
t95 = -mrSges(2,2) + mrSges(3,3);
t72 = cos(qJ(3));
t92 = qJD(1) * qJD(3);
t60 = t72 * t92;
t68 = sin(qJ(3));
t54 = -t68 * qJDD(1) - t60;
t87 = t68 * t92;
t55 = t72 * qJDD(1) - t87;
t75 = qJD(1) ^ 2;
t77 = (t97 * t75) - t99;
t23 = (-t55 + t87) * pkin(7) + (-t54 + t60) * pkin(3) + t77;
t53 = (pkin(3) * t68 - pkin(7) * t72) * qJD(1);
t62 = t68 * qJD(1);
t74 = qJD(3) ^ 2;
t89 = t69 * g(1) - t73 * g(2);
t81 = -t75 * qJ(2) + qJDD(2) - t89;
t40 = t97 * qJDD(1) + t81;
t88 = -t72 * g(3) + t68 * t40;
t26 = -t74 * pkin(3) + qJDD(3) * pkin(7) - t53 * t62 + t88;
t67 = sin(qJ(4));
t71 = cos(qJ(4));
t94 = t67 * t23 + t71 * t26;
t93 = qJD(1) * t72;
t59 = t62 + qJD(4);
t52 = (mrSges(4,1) * t68 + mrSges(4,2) * t72) * qJD(1);
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t93;
t50 = t71 * qJD(3) - t67 * t93;
t30 = t50 * qJD(4) + t67 * qJDD(3) + t71 * t55;
t49 = qJDD(4) - t54;
t51 = t67 * qJD(3) + t71 * t93;
t86 = t71 * t23 - t67 * t26;
t12 = (t50 * t59 - t30) * pkin(8) + (t50 * t51 + t49) * pkin(4) + t86;
t29 = -t51 * qJD(4) + t71 * qJDD(3) - t67 * t55;
t38 = t59 * pkin(4) - t51 * pkin(8);
t48 = t50 ^ 2;
t13 = -t48 * pkin(4) + t29 * pkin(8) - t59 * t38 + t94;
t66 = sin(qJ(5));
t70 = cos(qJ(5));
t31 = t70 * t50 - t66 * t51;
t18 = t31 * qJD(5) + t66 * t29 + t70 * t30;
t32 = t66 * t50 + t70 * t51;
t20 = -t31 * mrSges(6,1) + t32 * mrSges(6,2);
t58 = qJD(5) + t59;
t27 = -t58 * mrSges(6,2) + t31 * mrSges(6,3);
t44 = qJDD(5) + t49;
t10 = m(6) * (t70 * t12 - t66 * t13) - t18 * mrSges(6,3) + t44 * mrSges(6,1) - t32 * t20 + t58 * t27;
t17 = -t32 * qJD(5) + t70 * t29 - t66 * t30;
t28 = t58 * mrSges(6,1) - t32 * mrSges(6,3);
t11 = m(6) * (t66 * t12 + t70 * t13) + t17 * mrSges(6,3) - t44 * mrSges(6,2) + t31 * t20 - t58 * t28;
t33 = -t50 * mrSges(5,1) + t51 * mrSges(5,2);
t34 = -t59 * mrSges(5,2) + t50 * mrSges(5,3);
t7 = m(5) * t86 + t49 * mrSges(5,1) - t30 * mrSges(5,3) + t70 * t10 + t66 * t11 - t51 * t33 + t59 * t34;
t35 = t59 * mrSges(5,1) - t51 * mrSges(5,3);
t8 = m(5) * t94 - t49 * mrSges(5,2) + t29 * mrSges(5,3) - t66 * t10 + t70 * t11 + t50 * t33 - t59 * t35;
t4 = m(4) * t88 - qJDD(3) * mrSges(4,2) + t54 * mrSges(4,3) - qJD(3) * t57 - t52 * t62 - t67 * t7 + t71 * t8;
t56 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t62;
t84 = t68 * g(3) + t72 * t40;
t25 = -qJDD(3) * pkin(3) - t74 * pkin(7) + t53 * t93 - t84;
t80 = t17 * mrSges(6,1) + t31 * t27 - m(6) * (-t29 * pkin(4) - t48 * pkin(8) + t51 * t38 + t25) - t18 * mrSges(6,2) - t32 * t28;
t76 = m(5) * t25 - t29 * mrSges(5,1) + t30 * mrSges(5,2) - t50 * t34 + t51 * t35 - t80;
t9 = m(4) * t84 + qJDD(3) * mrSges(4,1) - t55 * mrSges(4,3) + qJD(3) * t56 - t52 * t93 - t76;
t90 = t72 * t4 - t68 * t9;
t82 = -m(3) * (-qJDD(1) * pkin(1) + t81) - t68 * t4 - t72 * t9;
t79 = m(4) * t77 - t54 * mrSges(4,1) + t55 * mrSges(4,2) + t56 * t62 + t57 * t93 + t67 * t8 + t71 * t7;
t78 = -m(3) * (t75 * pkin(1) + t99) + t79;
t2 = m(2) * t85 + t95 * qJDD(1) - (t96 * t75) + t78;
t1 = m(2) * t89 + t96 * qJDD(1) + t95 * t75 + t82;
t3 = [-m(1) * g(1) - t69 * t1 + t73 * t2, t2, -m(3) * g(3) + t90, t4, t8, t11; -m(1) * g(2) + t73 * t1 + t69 * t2, t1, -(t75 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t78, t9, t7, t10; (-m(1) + t98) * g(3) + t90, t98 * g(3) + t90, qJDD(1) * mrSges(3,2) - t75 * mrSges(3,3) - t82, t79, t76, -t80;];
f_new = t3;
