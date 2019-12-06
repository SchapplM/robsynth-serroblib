% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP6
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
%   pkin=[a2,a3,a4,a5,d2,d3,d4,theta1]';
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
% Datum: 2019-12-05 16:53
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP6_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:50:50
% EndTime: 2019-12-05 16:50:52
% DurationCPUTime: 0.67s
% Computational Cost: add. (5890->127), mult. (11602->165), div. (0->0), fcn. (7149->8), ass. (0->65)
t66 = sin(qJ(3));
t68 = cos(qJ(3));
t84 = qJD(2) * qJD(3);
t81 = t68 * t84;
t47 = t66 * qJDD(2) + t81;
t48 = t68 * qJDD(2) - t66 * t84;
t86 = qJD(2) * t66;
t51 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t86;
t85 = qJD(2) * t68;
t52 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t85;
t70 = qJD(2) ^ 2;
t65 = sin(qJ(4));
t93 = cos(qJ(4));
t39 = (t65 * t68 + t93 * t66) * qJD(2);
t23 = t39 * qJD(4) + t65 * t47 - t93 * t48;
t38 = t65 * t86 - t93 * t85;
t24 = -t38 * qJD(4) + t93 * t47 + t65 * t48;
t61 = qJD(3) + qJD(4);
t34 = -t61 * mrSges(5,2) - t38 * mrSges(5,3);
t35 = t61 * mrSges(5,1) - t39 * mrSges(5,3);
t53 = qJD(3) * pkin(3) - pkin(7) * t86;
t62 = t68 ^ 2;
t64 = sin(pkin(8));
t87 = cos(pkin(8));
t50 = -t87 * g(1) - t64 * g(2);
t63 = -g(3) + qJDD(1);
t67 = sin(qJ(2));
t69 = cos(qJ(2));
t79 = -t67 * t50 + t69 * t63;
t76 = -qJDD(2) * pkin(2) - t79;
t73 = -t48 * pkin(3) + t53 * t86 + (-pkin(7) * t62 - pkin(6)) * t70 + t76;
t36 = -t61 * mrSges(6,1) + t39 * mrSges(6,2);
t37 = -t38 * mrSges(6,2) + t61 * mrSges(6,3);
t74 = t24 * mrSges(6,3) + t39 * t36 - m(6) * (-0.2e1 * qJD(5) * t39 + (t38 * t61 - t24) * qJ(5) + (t39 * t61 + t23) * pkin(4) + t73) - t23 * mrSges(6,1) - t38 * t37;
t72 = m(5) * t73 + t23 * mrSges(5,1) + t24 * mrSges(5,2) + t38 * t34 + t39 * t35 - t74;
t95 = (t66 * t51 - t68 * t52) * qJD(2) + m(4) * (-t70 * pkin(6) + t76) - t48 * mrSges(4,1) + t47 * mrSges(4,2) + t72;
t27 = t38 * pkin(4) - t39 * qJ(5);
t59 = t61 ^ 2;
t60 = qJDD(3) + qJDD(4);
t88 = t69 * t50 + t67 * t63;
t33 = -t70 * pkin(2) + qJDD(2) * pkin(6) + t88;
t49 = t64 * g(1) - t87 * g(2);
t80 = -t66 * t33 - t68 * t49;
t17 = (-t47 + t81) * pkin(7) + (t66 * t68 * t70 + qJDD(3)) * pkin(3) + t80;
t89 = t68 * t33 - t66 * t49;
t18 = -t62 * t70 * pkin(3) + t48 * pkin(7) - qJD(3) * t53 + t89;
t75 = t93 * t17 - t65 * t18;
t94 = m(6) * (-t60 * pkin(4) - t59 * qJ(5) + t39 * t27 + qJDD(5) - t75);
t92 = -mrSges(5,3) - mrSges(6,2);
t91 = t65 * t17 + t93 * t18;
t28 = t38 * mrSges(6,1) - t39 * mrSges(6,3);
t90 = -t38 * mrSges(5,1) - t39 * mrSges(5,2) - t28;
t10 = m(5) * t75 - t94 + (t34 + t37) * t61 + (mrSges(5,1) + mrSges(6,1)) * t60 + t90 * t39 + t92 * t24;
t46 = (-mrSges(4,1) * t68 + mrSges(4,2) * t66) * qJD(2);
t82 = m(6) * (-t59 * pkin(4) + t60 * qJ(5) + 0.2e1 * qJD(5) * t61 - t38 * t27 + t91) + t61 * t36 + t60 * mrSges(6,3);
t9 = m(5) * t91 - t60 * mrSges(5,2) + t92 * t23 - t61 * t35 + t90 * t38 + t82;
t5 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t47 * mrSges(4,3) + qJD(3) * t52 + t93 * t10 - t46 * t86 + t65 * t9;
t6 = m(4) * t89 - qJDD(3) * mrSges(4,2) + t48 * mrSges(4,3) - qJD(3) * t51 - t65 * t10 + t46 * t85 + t93 * t9;
t3 = m(3) * t88 - t70 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t66 * t5 + t68 * t6;
t8 = m(3) * t79 + qJDD(2) * mrSges(3,1) - t70 * mrSges(3,2) - t95;
t83 = m(2) * t63 + t67 * t3 + t69 * t8;
t78 = -t68 * t5 - t66 * t6;
t4 = (m(2) + m(3)) * t49 + t78;
t1 = m(2) * t50 + t69 * t3 - t67 * t8;
t2 = [-m(1) * g(1) + t87 * t1 - t64 * t4, t1, t3, t6, t9, -t23 * mrSges(6,2) - t38 * t28 + t82; -m(1) * g(2) + t64 * t1 + t87 * t4, t4, t8, t5, t10, -t74; -m(1) * g(3) + t83, t83, -m(3) * t49 - t78, t95, t72, -t60 * mrSges(6,1) + t24 * mrSges(6,2) + t39 * t28 - t61 * t37 + t94;];
f_new = t2;
