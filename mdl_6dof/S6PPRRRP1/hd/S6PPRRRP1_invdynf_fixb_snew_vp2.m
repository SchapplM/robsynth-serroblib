% Calculate vector of cutting forces with Newton-Euler
% S6PPRRRP1
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
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-05-04 20:30
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRRP1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRP1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:26:59
% EndTime: 2019-05-04 20:27:04
% DurationCPUTime: 2.34s
% Computational Cost: add. (33848->142), mult. (60846->188), div. (0->0), fcn. (46507->14), ass. (0->82)
t71 = sin(pkin(11));
t75 = cos(pkin(11));
t64 = -g(1) * t75 - g(2) * t71;
t70 = sin(pkin(12));
t74 = cos(pkin(12));
t63 = g(1) * t71 - g(2) * t75;
t69 = -g(3) + qJDD(1);
t73 = sin(pkin(6));
t77 = cos(pkin(6));
t89 = t63 * t77 + t69 * t73;
t36 = -t70 * t64 + t89 * t74;
t52 = -t63 * t73 + t69 * t77 + qJDD(2);
t72 = sin(pkin(7));
t76 = cos(pkin(7));
t116 = t36 * t76 + t52 * t72;
t82 = cos(qJ(4));
t102 = qJD(3) * t82;
t85 = qJD(3) ^ 2;
t37 = t74 * t64 + t89 * t70;
t80 = sin(qJ(3));
t83 = cos(qJ(3));
t97 = t116 * t80 + t83 * t37;
t29 = -pkin(3) * t85 + qJDD(3) * pkin(9) + t97;
t31 = -t36 * t72 + t52 * t76;
t79 = sin(qJ(4));
t105 = t82 * t29 + t79 * t31;
t103 = qJD(3) * t79;
t78 = sin(qJ(5));
t81 = cos(qJ(5));
t57 = qJD(4) * t81 - t78 * t103;
t101 = qJD(3) * qJD(4);
t94 = t82 * t101;
t61 = qJDD(3) * t79 + t94;
t42 = qJD(5) * t57 + qJDD(4) * t78 + t61 * t81;
t67 = qJD(5) - t102;
t47 = -mrSges(7,2) * t67 + mrSges(7,3) * t57;
t95 = t79 * t101;
t62 = qJDD(3) * t82 - t95;
t56 = qJDD(5) - t62;
t58 = qJD(4) * t78 + t81 * t103;
t60 = (-pkin(4) * t82 - pkin(10) * t79) * qJD(3);
t84 = qJD(4) ^ 2;
t22 = -pkin(4) * t84 + qJDD(4) * pkin(10) + t60 * t102 + t105;
t114 = t116 * t83 - t80 * t37;
t28 = -qJDD(3) * pkin(3) - t85 * pkin(9) - t114;
t25 = (-t61 - t94) * pkin(10) + (-t62 + t95) * pkin(4) + t28;
t93 = -t78 * t22 + t81 * t25;
t100 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t67 - t42) * qJ(6) + (t57 * t58 + t56) * pkin(5) + t93) + t67 * t47 + t56 * mrSges(7,1);
t44 = -mrSges(7,1) * t57 + mrSges(7,2) * t58;
t45 = -mrSges(6,1) * t57 + mrSges(6,2) * t58;
t48 = -mrSges(6,2) * t67 + mrSges(6,3) * t57;
t13 = m(6) * t93 + t56 * mrSges(6,1) + t67 * t48 + (-t45 - t44) * t58 + (-mrSges(6,3) - mrSges(7,3)) * t42 + t100;
t106 = t81 * t22 + t78 * t25;
t41 = -qJD(5) * t58 + qJDD(4) * t81 - t61 * t78;
t50 = mrSges(7,1) * t67 - mrSges(7,3) * t58;
t51 = mrSges(6,1) * t67 - mrSges(6,3) * t58;
t49 = pkin(5) * t67 - qJ(6) * t58;
t55 = t57 ^ 2;
t99 = m(7) * (-pkin(5) * t55 + qJ(6) * t41 + 0.2e1 * qJD(6) * t57 - t49 * t67 + t106) + t57 * t44 + t41 * mrSges(7,3);
t14 = m(6) * t106 + t41 * mrSges(6,3) + t57 * t45 + (-t51 - t50) * t67 + (-mrSges(6,2) - mrSges(7,2)) * t56 + t99;
t59 = (-mrSges(5,1) * t82 + mrSges(5,2) * t79) * qJD(3);
t65 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t103;
t12 = m(5) * t105 - qJDD(4) * mrSges(5,2) + t62 * mrSges(5,3) - qJD(4) * t65 + t59 * t102 - t78 * t13 + t81 * t14;
t92 = -t79 * t29 + t31 * t82;
t21 = -qJDD(4) * pkin(4) - pkin(10) * t84 + t60 * t103 - t92;
t98 = m(7) * (-pkin(5) * t41 - qJ(6) * t55 + t49 * t58 + qJDD(6) + t21) + t42 * mrSges(7,2) + t58 * t50;
t112 = m(6) * t21 + t42 * mrSges(6,2) - (t48 + t47) * t57 - (mrSges(6,1) + mrSges(7,1)) * t41 + t58 * t51 + t98;
t66 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t102;
t15 = m(5) * t92 + qJDD(4) * mrSges(5,1) - t61 * mrSges(5,3) + qJD(4) * t66 - t59 * t103 - t112;
t10 = m(4) * t31 + t12 * t79 + t15 * t82;
t113 = m(5) * t28 - t62 * mrSges(5,1) + t61 * mrSges(5,2) + t81 * t13 + t78 * t14 + (t65 * t79 - t66 * t82) * qJD(3);
t11 = m(4) * t114 + qJDD(3) * mrSges(4,1) - t85 * mrSges(4,2) - t113;
t9 = m(4) * t97 - t85 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t82 * t12 - t79 * t15;
t91 = t11 * t83 + t80 * t9;
t4 = m(3) * t36 - t72 * t10 + t91 * t76;
t8 = m(3) * t37 - t11 * t80 + t83 * t9;
t115 = t4 * t74 + t70 * t8;
t6 = m(3) * t52 + t76 * t10 + t91 * t72;
t96 = m(2) * t69 + t115 * t73 + t77 * t6;
t2 = m(2) * t64 - t4 * t70 + t74 * t8;
t1 = m(2) * t63 + t115 * t77 - t73 * t6;
t3 = [-m(1) * g(1) - t1 * t71 + t2 * t75, t2, t8, t9, t12, t14, -t56 * mrSges(7,2) - t67 * t50 + t99; -m(1) * g(2) + t1 * t75 + t2 * t71, t1, t4, t11, t15, t13, -t42 * mrSges(7,3) - t58 * t44 + t100; -m(1) * g(3) + t96, t96, t6, t10, t113, t112, -t41 * mrSges(7,1) - t57 * t47 + t98;];
f_new  = t3;
