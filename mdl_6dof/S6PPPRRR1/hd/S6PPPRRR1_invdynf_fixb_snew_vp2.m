% Calculate vector of cutting forces with Newton-Euler
% S6PPPRRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% qJDD [6x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [14x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,alpha4,d4,d5,d6,theta1,theta2,theta3]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x7]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-04 19:36
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(14,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [14 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: pkin has to be [14x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 19:32:01
% EndTime: 2019-05-04 19:32:07
% DurationCPUTime: 5.80s
% Computational Cost: add. (97101->117), mult. (164579->170), div. (0->0), fcn. (146190->18), ass. (0->83)
t58 = sin(pkin(12));
t64 = cos(pkin(12));
t50 = -t64 * g(1) - t58 * g(2);
t57 = sin(pkin(13));
t63 = cos(pkin(13));
t49 = t58 * g(1) - t64 * g(2);
t55 = -g(3) + qJDD(1);
t61 = sin(pkin(6));
t67 = cos(pkin(6));
t79 = t49 * t67 + t55 * t61;
t35 = t63 * t50 + t79 * t57;
t56 = sin(pkin(14));
t62 = cos(pkin(14));
t34 = -t57 * t50 + t79 * t63;
t41 = -t61 * t49 + t67 * t55 + qJDD(2);
t60 = sin(pkin(7));
t66 = cos(pkin(7));
t80 = t34 * t66 + t41 * t60;
t30 = -t56 * t35 + t80 * t62;
t33 = -t60 * t34 + t66 * t41 + qJDD(3);
t59 = sin(pkin(8));
t65 = cos(pkin(8));
t100 = t30 * t65 + t33 * t59;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t46 = (-pkin(5) * t72 - pkin(11) * t69) * qJD(4);
t74 = qJD(5) ^ 2;
t89 = t72 * qJD(4);
t75 = qJD(4) ^ 2;
t31 = t62 * t35 + t80 * t56;
t70 = sin(qJ(4));
t73 = cos(qJ(4));
t87 = t100 * t70 + t73 * t31;
t24 = -t75 * pkin(4) + qJDD(4) * pkin(10) + t87;
t26 = -t59 * t30 + t65 * t33;
t91 = t72 * t24 + t69 * t26;
t20 = -t74 * pkin(5) + qJDD(5) * pkin(11) + t46 * t89 + t91;
t98 = t100 * t73 - t70 * t31;
t23 = -qJDD(4) * pkin(4) - t75 * pkin(10) - t98;
t88 = qJD(4) * qJD(5);
t84 = t72 * t88;
t47 = t69 * qJDD(4) + t84;
t85 = t69 * t88;
t48 = t72 * qJDD(4) - t85;
t21 = (-t47 - t84) * pkin(11) + (-t48 + t85) * pkin(5) + t23;
t68 = sin(qJ(6));
t71 = cos(qJ(6));
t90 = qJD(4) * t69;
t43 = t71 * qJD(5) - t68 * t90;
t37 = t43 * qJD(6) + t68 * qJDD(5) + t71 * t47;
t44 = t68 * qJD(5) + t71 * t90;
t38 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t53 = qJD(6) - t89;
t39 = -t53 * mrSges(7,2) + t43 * mrSges(7,3);
t42 = qJDD(6) - t48;
t17 = m(7) * (-t68 * t20 + t71 * t21) - t37 * mrSges(7,3) + t42 * mrSges(7,1) - t44 * t38 + t53 * t39;
t36 = -t44 * qJD(6) + t71 * qJDD(5) - t68 * t47;
t40 = t53 * mrSges(7,1) - t44 * mrSges(7,3);
t18 = m(7) * (t71 * t20 + t68 * t21) + t36 * mrSges(7,3) - t42 * mrSges(7,2) + t43 * t38 - t53 * t40;
t45 = (-mrSges(6,1) * t72 + mrSges(6,2) * t69) * qJD(4);
t51 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t90;
t15 = m(6) * t91 - qJDD(5) * mrSges(6,2) + t48 * mrSges(6,3) - qJD(5) * t51 - t68 * t17 + t71 * t18 + t45 * t89;
t52 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t89;
t92 = t72 * t26;
t76 = m(7) * (-qJDD(5) * pkin(5) - t74 * pkin(11) - t92 + (qJD(4) * t46 + t24) * t69) - t36 * mrSges(7,1) + t37 * mrSges(7,2) - t43 * t39 + t44 * t40;
t16 = m(6) * (-t69 * t24 + t92) - t47 * mrSges(6,3) + qJDD(5) * mrSges(6,1) - t45 * t90 + qJD(5) * t52 - t76;
t13 = m(5) * t26 + t69 * t15 + t72 * t16;
t12 = m(5) * t87 - t75 * mrSges(5,1) - qJDD(4) * mrSges(5,2) + t72 * t15 - t69 * t16;
t97 = m(6) * t23 - t48 * mrSges(6,1) + t47 * mrSges(6,2) + t71 * t17 + t68 * t18 + (t69 * t51 - t72 * t52) * qJD(4);
t14 = m(5) * t98 + qJDD(4) * mrSges(5,1) - t75 * mrSges(5,2) - t97;
t82 = t12 * t70 + t14 * t73;
t10 = m(4) * t33 + t65 * t13 + t82 * t59;
t11 = m(4) * t31 + t73 * t12 - t70 * t14;
t9 = m(4) * t30 - t59 * t13 + t82 * t65;
t83 = t11 * t56 + t62 * t9;
t4 = m(3) * t34 - t60 * t10 + t83 * t66;
t8 = m(3) * t35 + t62 * t11 - t56 * t9;
t99 = t4 * t63 + t57 * t8;
t6 = m(3) * t41 + t66 * t10 + t83 * t60;
t86 = m(2) * t55 + t67 * t6 + t99 * t61;
t2 = m(2) * t50 - t57 * t4 + t63 * t8;
t1 = m(2) * t49 - t61 * t6 + t99 * t67;
t3 = [-m(1) * g(1) - t58 * t1 + t64 * t2, t2, t8, t11, t12, t15, t18; -m(1) * g(2) + t64 * t1 + t58 * t2, t1, t4, t9, t14, t16, t17; -m(1) * g(3) + t86, t86, t6, t10, t13, t97, t76;];
f_new  = t3;
