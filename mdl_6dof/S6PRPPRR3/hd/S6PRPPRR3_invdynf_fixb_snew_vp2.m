% Calculate vector of cutting forces with Newton-Euler
% S6PRPPRR3
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-05-04 22:02
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 21:59:57
% EndTime: 2019-05-04 21:59:59
% DurationCPUTime: 1.05s
% Computational Cost: add. (13821->125), mult. (23929->165), div. (0->0), fcn. (14298->12), ass. (0->77)
t61 = sin(pkin(10));
t64 = cos(pkin(10));
t49 = -g(1) * t64 - g(2) * t61;
t58 = -g(3) + qJDD(1);
t62 = sin(pkin(6));
t68 = sin(qJ(2));
t71 = cos(qJ(2));
t48 = g(1) * t61 - g(2) * t64;
t65 = cos(pkin(6));
t96 = t48 * t65;
t101 = -t68 * t49 + (t58 * t62 + t96) * t71;
t67 = sin(qJ(5));
t70 = cos(qJ(5));
t45 = (pkin(5) * t70 + pkin(9) * t67) * qJD(2);
t72 = qJD(5) ^ 2;
t90 = qJD(2) * t70;
t73 = qJD(2) ^ 2;
t95 = t62 * t68;
t88 = t49 * t71 + t58 * t95 + t68 * t96;
t83 = qJDD(2) * qJ(3) + 0.2e1 * qJD(3) * qJD(2) + t88;
t99 = -pkin(2) - pkin(3);
t25 = t73 * t99 + t83;
t76 = -qJ(3) * t73 + qJDD(3) - t101;
t27 = qJDD(2) * t99 + t76;
t60 = sin(pkin(11));
t63 = cos(pkin(11));
t92 = t25 * t63 + t27 * t60;
t21 = -pkin(4) * t73 - qJDD(2) * pkin(8) + t92;
t37 = -t48 * t62 + t58 * t65;
t35 = qJDD(4) - t37;
t93 = t21 * t70 + t35 * t67;
t17 = -pkin(5) * t72 + qJDD(5) * pkin(9) - t45 * t90 + t93;
t84 = -t25 * t60 + t27 * t63;
t20 = qJDD(2) * pkin(4) - pkin(8) * t73 - t84;
t89 = qJD(2) * qJD(5);
t85 = t70 * t89;
t46 = -qJDD(2) * t67 - t85;
t86 = t67 * t89;
t47 = -qJDD(2) * t70 + t86;
t18 = (-t46 + t85) * pkin(9) + (-t47 - t86) * pkin(5) + t20;
t66 = sin(qJ(6));
t69 = cos(qJ(6));
t91 = qJD(2) * t67;
t42 = qJD(5) * t69 + t66 * t91;
t30 = qJD(6) * t42 + qJDD(5) * t66 + t46 * t69;
t43 = qJD(5) * t66 - t69 * t91;
t31 = -mrSges(7,1) * t42 + mrSges(7,2) * t43;
t53 = qJD(6) + t90;
t33 = -mrSges(7,2) * t53 + mrSges(7,3) * t42;
t39 = qJDD(6) - t47;
t14 = m(7) * (-t17 * t66 + t18 * t69) - t30 * mrSges(7,3) + t39 * mrSges(7,1) - t43 * t31 + t53 * t33;
t29 = -qJD(6) * t43 + qJDD(5) * t69 - t46 * t66;
t34 = mrSges(7,1) * t53 - mrSges(7,3) * t43;
t15 = m(7) * (t17 * t69 + t18 * t66) + t29 * mrSges(7,3) - t39 * mrSges(7,2) + t42 * t31 - t53 * t34;
t50 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t91;
t51 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t90;
t100 = m(6) * t20 - t47 * mrSges(6,1) + t46 * mrSges(6,2) - (t50 * t67 - t51 * t70) * qJD(2) + t69 * t14 + t66 * t15;
t11 = m(5) * t84 - qJDD(2) * mrSges(5,1) - t73 * mrSges(5,2) - t100;
t44 = (mrSges(6,1) * t70 - mrSges(6,2) * t67) * qJD(2);
t12 = m(6) * t93 - qJDD(5) * mrSges(6,2) + mrSges(6,3) * t47 - qJD(5) * t50 - t14 * t66 + t15 * t69 - t44 * t90;
t97 = t35 * t70;
t74 = m(7) * (-qJDD(5) * pkin(5) - pkin(9) * t72 - t97 + (-qJD(2) * t45 + t21) * t67) - t29 * mrSges(7,1) + t30 * mrSges(7,2) - t42 * t33 + t43 * t34;
t13 = m(6) * (-t21 * t67 + t97) - t46 * mrSges(6,3) + qJDD(5) * mrSges(6,1) + t44 * t91 + qJD(5) * t51 - t74;
t8 = m(5) * t92 - mrSges(5,1) * t73 + qJDD(2) * mrSges(5,2) + t12 * t70 - t13 * t67;
t79 = -m(4) * (-qJDD(2) * pkin(2) + t76) - t63 * t11 - t60 * t8;
t94 = mrSges(3,1) + mrSges(4,1);
t6 = m(3) * t101 + (-mrSges(3,2) + mrSges(4,3)) * t73 + t94 * qJDD(2) + t79;
t98 = t6 * t71;
t78 = m(5) * t35 + t12 * t67 + t13 * t70;
t77 = m(4) * t37 - t78;
t10 = m(3) * t37 + t77;
t80 = -t60 * t11 + t63 * t8 + m(4) * (-pkin(2) * t73 + t83) + qJDD(2) * mrSges(4,3);
t5 = m(3) * t88 - qJDD(2) * mrSges(3,2) - t73 * t94 + t80;
t87 = m(2) * t58 + t10 * t65 + t5 * t95 + t62 * t98;
t2 = m(2) * t49 + t5 * t71 - t6 * t68;
t1 = m(2) * t48 - t10 * t62 + (t5 * t68 + t98) * t65;
t3 = [-m(1) * g(1) - t1 * t61 + t2 * t64, t2, t5, -t73 * mrSges(4,1) + t80, t8, t12, t15; -m(1) * g(2) + t1 * t64 + t2 * t61, t1, t6, t77, t11, t13, t14; -m(1) * g(3) + t87, t87, t10, -qJDD(2) * mrSges(4,1) - t73 * mrSges(4,3) - t79, t78, t100, t74;];
f_new  = t3;
