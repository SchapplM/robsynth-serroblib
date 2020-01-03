% Calculate vector of cutting forces with Newton-Euler
% S5RPRPP4
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta4]';
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
% Datum: 2019-12-31 18:15
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPP4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP4_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:14:19
% EndTime: 2019-12-31 18:14:20
% DurationCPUTime: 0.57s
% Computational Cost: add. (3844->133), mult. (8123->166), div. (0->0), fcn. (4493->6), ass. (0->63)
t63 = sin(qJ(1));
t65 = cos(qJ(1));
t78 = -t65 * g(1) - t63 * g(2);
t75 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t78;
t97 = -2 * qJD(4);
t60 = sin(pkin(7));
t61 = cos(pkin(7));
t62 = sin(qJ(3));
t64 = cos(qJ(3));
t40 = (t60 * t64 + t61 * t62) * qJD(1);
t87 = qJD(1) * t64;
t88 = qJD(1) * t62;
t41 = -t60 * t88 + t61 * t87;
t20 = t40 * pkin(4) - t41 * qJ(5);
t66 = qJD(3) ^ 2;
t86 = qJD(1) * qJD(3);
t79 = t62 * t86;
t49 = t64 * qJDD(1) - t79;
t67 = qJD(1) ^ 2;
t81 = t63 * g(1) - t65 * g(2);
t73 = -t67 * qJ(2) + qJDD(2) - t81;
t94 = -pkin(1) - pkin(6);
t36 = t94 * qJDD(1) + t73;
t89 = t62 * g(3) + t64 * t36;
t14 = (-t49 - t79) * qJ(4) + (-t62 * t64 * t67 + qJDD(3)) * pkin(3) + t89;
t48 = -t62 * qJDD(1) - t64 * t86;
t51 = qJD(3) * pkin(3) - qJ(4) * t87;
t59 = t62 ^ 2;
t80 = -t64 * g(3) + t62 * t36;
t15 = -t59 * t67 * pkin(3) + t48 * qJ(4) - qJD(3) * t51 + t80;
t77 = t61 * t14 - t60 * t15;
t96 = m(6) * (-qJDD(3) * pkin(4) - t66 * qJ(5) + qJDD(5) + ((2 * qJD(4)) + t20) * t41 - t77);
t95 = -m(2) - m(3);
t93 = (mrSges(2,1) - mrSges(3,2));
t92 = -mrSges(2,2) + mrSges(3,3);
t91 = -mrSges(5,3) - mrSges(6,2);
t21 = t40 * mrSges(6,1) - t41 * mrSges(6,3);
t90 = -t40 * mrSges(5,1) - t41 * mrSges(5,2) - t21;
t34 = -qJD(3) * mrSges(6,1) + t41 * mrSges(6,2);
t83 = t60 * t14 + t61 * t15 + t40 * t97;
t84 = qJD(3) * t34 + qJDD(3) * mrSges(6,3) + m(6) * (-t66 * pkin(4) + qJDD(3) * qJ(5) + (2 * qJD(5) * qJD(3)) - t40 * t20 + t83);
t47 = (mrSges(4,1) * t62 + mrSges(4,2) * t64) * qJD(1);
t50 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t88;
t25 = -t61 * t48 + t60 * t49;
t33 = qJD(3) * mrSges(5,1) - t41 * mrSges(5,3);
t6 = m(5) * t83 - qJDD(3) * mrSges(5,2) - qJD(3) * t33 + t91 * t25 + t90 * t40 + t84;
t26 = t60 * t48 + t61 * t49;
t32 = -qJD(3) * mrSges(5,2) - t40 * mrSges(5,3);
t35 = -t40 * mrSges(6,2) + qJD(3) * mrSges(6,3);
t7 = m(5) * t77 - t96 + (m(5) * t97 + t90) * t41 + t91 * t26 + (mrSges(5,1) + mrSges(6,1)) * qJDD(3) + (t32 + t35) * qJD(3);
t3 = m(4) * t89 + qJDD(3) * mrSges(4,1) - t49 * mrSges(4,3) + qJD(3) * t50 - t47 * t87 + t60 * t6 + t61 * t7;
t52 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t87;
t4 = m(4) * t80 - qJDD(3) * mrSges(4,2) + t48 * mrSges(4,3) - qJD(3) * t52 - t47 * t88 + t61 * t6 - t60 * t7;
t82 = -t62 * t3 + t64 * t4;
t74 = -m(3) * (-qJDD(1) * pkin(1) + t73) - t64 * t3 - t62 * t4;
t70 = -t48 * pkin(3) + qJDD(4) + t51 * t87 + (-qJ(4) * t59 + t94) * t67 + t75;
t72 = -t26 * mrSges(6,3) - t41 * t34 + m(6) * (-0.2e1 * qJD(5) * t41 + (qJD(3) * t40 - t26) * qJ(5) + (qJD(3) * t41 + t25) * pkin(4) + t70) + t40 * t35 + t25 * mrSges(6,1);
t71 = m(5) * t70 + t25 * mrSges(5,1) + t26 * mrSges(5,2) + t40 * t32 + t41 * t33 + t72;
t69 = -t48 * mrSges(4,1) + m(4) * (t94 * t67 + t75) + t50 * t88 + t52 * t87 + t49 * mrSges(4,2) + t71;
t68 = -m(3) * (t67 * pkin(1) - t75) + t69;
t5 = m(2) * t78 + t92 * qJDD(1) - (t93 * t67) + t68;
t1 = m(2) * t81 + t93 * qJDD(1) + t92 * t67 + t74;
t2 = [-m(1) * g(1) - t63 * t1 + t65 * t5, t5, -m(3) * g(3) + t82, t4, t6, -t25 * mrSges(6,2) - t40 * t21 + t84; -m(1) * g(2) + t65 * t1 + t63 * t5, t1, -(t67 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t68, t3, t7, t72; (-m(1) + t95) * g(3) + t82, t95 * g(3) + t82, qJDD(1) * mrSges(3,2) - t67 * mrSges(3,3) - t74, t69, t71, -qJDD(3) * mrSges(6,1) + t26 * mrSges(6,2) - qJD(3) * t35 + t41 * t21 + t96;];
f_new = t2;
