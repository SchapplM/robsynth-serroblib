% Calculate vector of cutting forces with Newton-Euler
% S5RPPRR1
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
%   pkin=[a2,a3,a4,a5,d1,d4,d5]';
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
% Datum: 2019-12-05 17:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR1_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:38:01
% EndTime: 2019-12-05 17:38:02
% DurationCPUTime: 0.46s
% Computational Cost: add. (2757->112), mult. (5263->137), div. (0->0), fcn. (2561->6), ass. (0->58)
t93 = 2 * qJD(1);
t61 = sin(qJ(1));
t64 = cos(qJ(1));
t74 = -t64 * g(1) - t61 * g(2);
t92 = qJDD(1) * qJ(2) + (qJD(2) * t93) + t74;
t65 = qJD(1) ^ 2;
t87 = t61 * g(1) - t64 * g(2);
t33 = -qJDD(1) * pkin(1) - t65 * qJ(2) + qJDD(2) - t87;
t70 = qJDD(1) * qJ(3) + (qJD(3) * t93) - t33;
t91 = -m(3) - m(4);
t90 = (mrSges(3,2) - mrSges(4,3));
t89 = (-mrSges(4,2) - mrSges(3,3));
t67 = qJDD(3) + (-pkin(1) - qJ(3)) * t65 + t92;
t25 = -qJDD(1) * pkin(6) + t67;
t60 = sin(qJ(4));
t63 = cos(qJ(4));
t88 = t60 * g(3) + t63 * t25;
t86 = qJD(1) * t60;
t85 = qJD(1) * t63;
t84 = -m(2) + t91;
t83 = qJD(1) * qJD(4);
t80 = mrSges(2,1) - t90;
t40 = (mrSges(5,1) * t60 + mrSges(5,2) * t63) * qJD(1);
t77 = t60 * t83;
t42 = t63 * qJDD(1) - t77;
t43 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t86;
t59 = sin(qJ(5));
t62 = cos(qJ(5));
t10 = (-t42 - t77) * pkin(7) + (-t60 * t63 * t65 + qJDD(4)) * pkin(4) + t88;
t41 = -t60 * qJDD(1) - t63 * t83;
t45 = qJD(4) * pkin(4) - pkin(7) * t85;
t57 = t60 ^ 2;
t78 = -t63 * g(3) + t60 * t25;
t11 = -t57 * t65 * pkin(4) + t41 * pkin(7) - qJD(4) * t45 + t78;
t34 = (-t59 * t63 - t60 * t62) * qJD(1);
t17 = t34 * qJD(5) + t59 * t41 + t62 * t42;
t35 = (-t59 * t60 + t62 * t63) * qJD(1);
t20 = -t34 * mrSges(6,1) + t35 * mrSges(6,2);
t52 = qJD(4) + qJD(5);
t30 = -t52 * mrSges(6,2) + t34 * mrSges(6,3);
t51 = qJDD(4) + qJDD(5);
t8 = m(6) * (t62 * t10 - t59 * t11) - t17 * mrSges(6,3) + t51 * mrSges(6,1) - t35 * t20 + t52 * t30;
t16 = -t35 * qJD(5) + t62 * t41 - t59 * t42;
t31 = t52 * mrSges(6,1) - t35 * mrSges(6,3);
t9 = m(6) * (t59 * t10 + t62 * t11) + t16 * mrSges(6,3) - t51 * mrSges(6,2) + t34 * t20 - t52 * t31;
t5 = m(5) * t88 + qJDD(4) * mrSges(5,1) - t42 * mrSges(5,3) + qJD(4) * t43 - t40 * t85 + t59 * t9 + t62 * t8;
t44 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t85;
t6 = m(5) * t78 - qJDD(4) * mrSges(5,2) + t41 * mrSges(5,3) - qJD(4) * t44 - t40 * t86 - t59 * t8 + t62 * t9;
t79 = -t60 * t5 + t63 * t6;
t76 = m(4) * t67 + qJDD(1) * mrSges(4,2) + t63 * t5 + t60 * t6;
t75 = m(6) * (t45 * t85 - t41 * pkin(4) + (-pkin(7) * t57 - pkin(6)) * t65 + t70) + t17 * mrSges(6,2) - t16 * mrSges(6,1) + t35 * t31 - t34 * t30;
t72 = m(3) * (t65 * pkin(1) - t92) - t76;
t69 = m(5) * (-t65 * pkin(6) + t70) + t42 * mrSges(5,2) - t41 * mrSges(5,1) + t44 * t85 + t43 * t86 + t75;
t68 = -m(4) * t70 - t69;
t66 = m(3) * t33 + t68;
t7 = m(2) * t87 - t66 + t80 * qJDD(1) + ((-mrSges(2,2) - t89) * t65);
t1 = m(2) * t74 + (-mrSges(2,2) + mrSges(3,3)) * qJDD(1) - t80 * t65 - t72;
t2 = [-m(1) * g(1) + t64 * t1 - t61 * t7, t1, t91 * g(3) + t79, -m(4) * g(3) + t79, t6, t9; -m(1) * g(2) + t61 * t1 + t64 * t7, t7, -qJDD(1) * mrSges(3,3) - (t90 * t65) + t72, -t65 * mrSges(4,2) - qJDD(1) * mrSges(4,3) + t68, t5, t8; (-m(1) + t84) * g(3) + t79, t84 * g(3) + t79, t90 * qJDD(1) + (t89 * t65) + t66, -t65 * mrSges(4,3) + t76, t69, t75;];
f_new = t2;
