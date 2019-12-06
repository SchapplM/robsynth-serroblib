% Calculate vector of cutting forces with Newton-Euler
% S5PRRRP5
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
% Datum: 2019-12-05 16:49
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5PRRRP5_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRP5_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:47:49
% EndTime: 2019-12-05 16:47:50
% DurationCPUTime: 0.72s
% Computational Cost: add. (6076->128), mult. (12235->166), div. (0->0), fcn. (7635->8), ass. (0->64)
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t87 = qJD(2) * qJD(3);
t82 = t70 * t87;
t50 = t67 * qJDD(2) + t82;
t51 = t70 * qJDD(2) - t67 * t87;
t89 = qJD(2) * t67;
t54 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t89;
t88 = qJD(2) * t70;
t55 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t88;
t72 = qJD(2) ^ 2;
t66 = sin(qJ(4));
t69 = cos(qJ(4));
t43 = (t66 * t70 + t67 * t69) * qJD(2);
t26 = -t43 * qJD(4) - t66 * t50 + t69 * t51;
t42 = (-t66 * t67 + t69 * t70) * qJD(2);
t27 = t42 * qJD(4) + t69 * t50 + t66 * t51;
t62 = qJD(3) + qJD(4);
t36 = -t62 * mrSges(6,2) + t42 * mrSges(6,3);
t37 = -t62 * mrSges(5,2) + t42 * mrSges(5,3);
t40 = t62 * mrSges(5,1) - t43 * mrSges(5,3);
t56 = qJD(3) * pkin(3) - pkin(7) * t89;
t63 = t70 ^ 2;
t65 = sin(pkin(8));
t90 = cos(pkin(8));
t53 = -t90 * g(1) - t65 * g(2);
t64 = -g(3) + qJDD(1);
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t79 = -t68 * t53 + t71 * t64;
t75 = -qJDD(2) * pkin(2) - t79;
t73 = -t51 * pkin(3) + t56 * t89 + (-pkin(7) * t63 - pkin(6)) * t72 + t75;
t38 = t62 * pkin(4) - t43 * qJ(5);
t39 = t62 * mrSges(6,1) - t43 * mrSges(6,3);
t41 = t42 ^ 2;
t83 = m(6) * (-t26 * pkin(4) - t41 * qJ(5) + t43 * t38 + qJDD(5) + t73) + t27 * mrSges(6,2) + t43 * t39;
t74 = m(5) * t73 + t27 * mrSges(5,2) - (mrSges(5,1) + mrSges(6,1)) * t26 + t43 * t40 - (t37 + t36) * t42 + t83;
t100 = t74 + (t67 * t54 - t70 * t55) * qJD(2) - t51 * mrSges(4,1) + t50 * mrSges(4,2) + m(4) * (-t72 * pkin(6) + t75);
t91 = t71 * t53 + t68 * t64;
t35 = -t72 * pkin(2) + qJDD(2) * pkin(6) + t91;
t52 = t65 * g(1) - t90 * g(2);
t80 = -t67 * t35 - t70 * t52;
t18 = (-t50 + t82) * pkin(7) + (t67 * t70 * t72 + qJDD(3)) * pkin(3) + t80;
t93 = t70 * t35 - t67 * t52;
t19 = -t63 * t72 * pkin(3) + t51 * pkin(7) - qJD(3) * t56 + t93;
t94 = t66 * t18 + t69 * t19;
t31 = -t42 * mrSges(5,1) + t43 * mrSges(5,2);
t61 = qJDD(3) + qJDD(4);
t30 = -t42 * mrSges(6,1) + t43 * mrSges(6,2);
t84 = m(6) * (-t41 * pkin(4) + t26 * qJ(5) + 0.2e1 * qJD(5) * t42 - t62 * t38 + t94) + t42 * t30 + t26 * mrSges(6,3);
t10 = m(5) * t94 + t26 * mrSges(5,3) + t42 * t31 + (-t40 - t39) * t62 + (-mrSges(5,2) - mrSges(6,2)) * t61 + t84;
t49 = (-mrSges(4,1) * t70 + mrSges(4,2) * t67) * qJD(2);
t81 = t69 * t18 - t66 * t19;
t85 = m(6) * (-0.2e1 * qJD(5) * t43 + (t42 * t62 - t27) * qJ(5) + (t42 * t43 + t61) * pkin(4) + t81) + t62 * t36 + t61 * mrSges(6,1);
t7 = m(5) * t81 + t61 * mrSges(5,1) + t62 * t37 + (-t31 - t30) * t43 + (-mrSges(5,3) - mrSges(6,3)) * t27 + t85;
t5 = m(4) * t80 + qJDD(3) * mrSges(4,1) - t50 * mrSges(4,3) + qJD(3) * t55 + t66 * t10 - t49 * t89 + t69 * t7;
t6 = m(4) * t93 - qJDD(3) * mrSges(4,2) + t51 * mrSges(4,3) - qJD(3) * t54 + t69 * t10 + t49 * t88 - t66 * t7;
t3 = m(3) * t91 - t72 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t67 * t5 + t70 * t6;
t9 = m(3) * t79 + qJDD(2) * mrSges(3,1) - t72 * mrSges(3,2) - t100;
t86 = m(2) * t64 + t68 * t3 + t71 * t9;
t78 = -t70 * t5 - t67 * t6;
t4 = (m(2) + m(3)) * t52 + t78;
t1 = m(2) * t53 + t71 * t3 - t68 * t9;
t2 = [-m(1) * g(1) + t90 * t1 - t65 * t4, t1, t3, t6, t10, -t61 * mrSges(6,2) - t62 * t39 + t84; -m(1) * g(2) + t65 * t1 + t90 * t4, t4, t9, t5, t7, -t27 * mrSges(6,3) - t43 * t30 + t85; -m(1) * g(3) + t86, t86, -m(3) * t52 - t78, t100, t74, -t26 * mrSges(6,1) - t42 * t36 + t83;];
f_new = t2;
