% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR15
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:37
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR15_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR15_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:36:40
% EndTime: 2019-12-31 18:36:42
% DurationCPUTime: 0.81s
% Computational Cost: add. (7994->137), mult. (16560->176), div. (0->0), fcn. (9821->8), ass. (0->73)
t67 = sin(qJ(1));
t70 = cos(qJ(1));
t82 = -t70 * g(1) - t67 * g(2);
t98 = -qJDD(1) * qJ(2) - (2 * qJD(2) * qJD(1)) - t82;
t97 = -m(2) - m(3);
t96 = (-pkin(1) - pkin(6));
t95 = (mrSges(2,1) - mrSges(3,2));
t94 = -mrSges(2,2) + mrSges(3,3);
t69 = cos(qJ(3));
t93 = qJD(1) * t69;
t66 = sin(qJ(3));
t92 = t66 * qJD(1);
t91 = qJD(1) * qJD(3);
t84 = t69 * t91;
t55 = t66 * qJDD(1) + t84;
t85 = t66 * t91;
t56 = t69 * qJDD(1) - t85;
t72 = qJD(1) ^ 2;
t74 = (t96 * t72) - t98;
t23 = (-t56 + t85) * qJ(4) + (t55 + t84) * pkin(3) + t74;
t53 = (pkin(3) * t66 - qJ(4) * t69) * qJD(1);
t71 = qJD(3) ^ 2;
t87 = t67 * g(1) - t70 * g(2);
t78 = -t72 * qJ(2) + qJDD(2) - t87;
t40 = t96 * qJDD(1) + t78;
t86 = -t69 * g(3) + t66 * t40;
t26 = -t71 * pkin(3) + qJDD(3) * qJ(4) - t53 * t92 + t86;
t63 = sin(pkin(8));
t64 = cos(pkin(8));
t49 = t64 * qJD(3) - t63 * t93;
t89 = 0.2e1 * qJD(4) * t49 + t63 * t23 + t64 * t26;
t54 = (mrSges(4,1) * t66 + mrSges(4,2) * t69) * qJD(1);
t58 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t93;
t38 = t63 * qJDD(3) + t64 * t56;
t50 = t63 * qJD(3) + t64 * t93;
t83 = -0.2e1 * qJD(4) * t50 + t64 * t23 - t63 * t26;
t12 = (t49 * t92 - t38) * pkin(7) + (t49 * t50 + t55) * pkin(4) + t83;
t37 = t64 * qJDD(3) - t63 * t56;
t39 = pkin(4) * t92 - t50 * pkin(7);
t48 = t49 ^ 2;
t13 = -t48 * pkin(4) + t37 * pkin(7) - t39 * t92 + t89;
t65 = sin(qJ(5));
t68 = cos(qJ(5));
t29 = t68 * t49 - t65 * t50;
t18 = t29 * qJD(5) + t65 * t37 + t68 * t38;
t30 = t65 * t49 + t68 * t50;
t20 = -t29 * mrSges(6,1) + t30 * mrSges(6,2);
t59 = qJD(5) + t92;
t27 = -t59 * mrSges(6,2) + t29 * mrSges(6,3);
t52 = qJDD(5) + t55;
t10 = m(6) * (t68 * t12 - t65 * t13) - t18 * mrSges(6,3) + t52 * mrSges(6,1) - t30 * t20 + t59 * t27;
t17 = -t30 * qJD(5) + t68 * t37 - t65 * t38;
t28 = t59 * mrSges(6,1) - t30 * mrSges(6,3);
t11 = m(6) * (t65 * t12 + t68 * t13) + t17 * mrSges(6,3) - t52 * mrSges(6,2) + t29 * t20 - t59 * t28;
t31 = -t49 * mrSges(5,1) + t50 * mrSges(5,2);
t35 = -mrSges(5,2) * t92 + t49 * mrSges(5,3);
t7 = m(5) * t83 + t55 * mrSges(5,1) - t38 * mrSges(5,3) + t68 * t10 + t65 * t11 - t50 * t31 + t35 * t92;
t36 = mrSges(5,1) * t92 - t50 * mrSges(5,3);
t8 = m(5) * t89 - t55 * mrSges(5,2) + t37 * mrSges(5,3) - t65 * t10 + t68 * t11 + t49 * t31 - t36 * t92;
t4 = m(4) * t86 - qJDD(3) * mrSges(4,2) - t55 * mrSges(4,3) - qJD(3) * t58 - t54 * t92 - t63 * t7 + t64 * t8;
t57 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t92;
t81 = t66 * g(3) + t69 * t40;
t25 = -qJDD(3) * pkin(3) - t71 * qJ(4) + t53 * t93 + qJDD(4) - t81;
t77 = t17 * mrSges(6,1) + t29 * t27 - m(6) * (-t37 * pkin(4) - t48 * pkin(7) + t50 * t39 + t25) - t18 * mrSges(6,2) - t30 * t28;
t73 = m(5) * t25 - t37 * mrSges(5,1) + t38 * mrSges(5,2) - t49 * t35 + t50 * t36 - t77;
t9 = m(4) * t81 + qJDD(3) * mrSges(4,1) - t56 * mrSges(4,3) + qJD(3) * t57 - t54 * t93 - t73;
t88 = t69 * t4 - t66 * t9;
t79 = -m(3) * (-qJDD(1) * pkin(1) + t78) - t66 * t4 - t69 * t9;
t76 = m(4) * t74 + t55 * mrSges(4,1) + t56 * mrSges(4,2) + t57 * t92 + t58 * t93 + t63 * t8 + t64 * t7;
t75 = -m(3) * (t72 * pkin(1) + t98) + t76;
t2 = m(2) * t82 + t94 * qJDD(1) - (t95 * t72) + t75;
t1 = m(2) * t87 + t95 * qJDD(1) + t94 * t72 + t79;
t3 = [-m(1) * g(1) - t67 * t1 + t70 * t2, t2, -m(3) * g(3) + t88, t4, t8, t11; -m(1) * g(2) + t70 * t1 + t67 * t2, t1, -(t72 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t75, t9, t7, t10; (-m(1) + t97) * g(3) + t88, t97 * g(3) + t88, qJDD(1) * mrSges(3,2) - t72 * mrSges(3,3) - t79, t76, t73, -t77;];
f_new = t3;
