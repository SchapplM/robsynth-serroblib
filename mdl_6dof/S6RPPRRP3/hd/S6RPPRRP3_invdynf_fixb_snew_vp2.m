% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRP3
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,theta2]';
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
% Datum: 2019-05-05 14:52
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRP3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRP3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:50:44
% EndTime: 2019-05-05 14:50:46
% DurationCPUTime: 0.75s
% Computational Cost: add. (7532->147), mult. (13685->177), div. (0->0), fcn. (7384->8), ass. (0->75)
t74 = sin(qJ(1));
t76 = cos(qJ(1));
t93 = t74 * g(1) - t76 * g(2);
t50 = qJDD(1) * pkin(1) + t93;
t78 = qJD(1) ^ 2;
t87 = -t76 * g(1) - t74 * g(2);
t52 = -t78 * pkin(1) + t87;
t70 = sin(pkin(9));
t71 = cos(pkin(9));
t102 = t70 * t50 + t71 * t52;
t114 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t102;
t75 = cos(qJ(4));
t100 = qJD(1) * t75;
t73 = sin(qJ(4));
t53 = (pkin(4) * t73 - pkin(8) * t75) * qJD(1);
t77 = qJD(4) ^ 2;
t112 = -pkin(2) - pkin(7);
t89 = t71 * t50 - t70 * t52;
t82 = -t78 * qJ(3) + qJDD(3) - t89;
t24 = qJDD(1) * t112 + t82;
t69 = -g(3) + qJDD(2);
t90 = t75 * t24 - t73 * t69;
t19 = -qJDD(4) * pkin(4) - t77 * pkin(8) + t53 * t100 - t90;
t110 = cos(qJ(5));
t72 = sin(qJ(5));
t49 = t72 * qJD(4) + t100 * t110;
t98 = qJD(1) * qJD(4);
t92 = t73 * t98;
t55 = t75 * qJDD(1) - t92;
t29 = t49 * qJD(5) - qJDD(4) * t110 + t72 * t55;
t48 = -qJD(4) * t110 + t100 * t72;
t30 = -t48 * qJD(5) + t72 * qJDD(4) + t110 * t55;
t99 = t73 * qJD(1);
t59 = qJD(5) + t99;
t35 = -t59 * mrSges(6,2) - t48 * mrSges(6,3);
t36 = t59 * mrSges(6,1) - t49 * mrSges(6,3);
t37 = -t59 * mrSges(7,1) + t49 * mrSges(7,2);
t38 = -t48 * mrSges(7,2) + t59 * mrSges(7,3);
t95 = m(7) * (-0.2e1 * qJD(6) * t49 + (t48 * t59 - t30) * qJ(6) + (t49 * t59 + t29) * pkin(5) + t19) + t29 * mrSges(7,1) + t48 * t38;
t113 = m(6) * t19 + t29 * mrSges(6,1) + (t36 - t37) * t49 + (mrSges(6,2) - mrSges(7,3)) * t30 + t48 * t35 + t95;
t32 = t48 * pkin(5) - t49 * qJ(6);
t91 = t75 * t98;
t54 = -t73 * qJDD(1) - t91;
t47 = qJDD(5) - t54;
t58 = t59 ^ 2;
t83 = t112 * t78 - t114;
t17 = (-t55 + t92) * pkin(8) + (-t54 + t91) * pkin(4) + t83;
t101 = t73 * t24 + t75 * t69;
t20 = -t77 * pkin(4) + qJDD(4) * pkin(8) - t53 * t99 + t101;
t85 = t110 * t17 - t72 * t20;
t111 = m(7) * (-t47 * pkin(5) - t58 * qJ(6) + t49 * t32 + qJDD(6) - t85);
t109 = mrSges(3,1) - mrSges(4,2);
t108 = -mrSges(3,2) + mrSges(4,3);
t106 = -mrSges(6,3) - mrSges(7,2);
t105 = t110 * t20 + t72 * t17;
t33 = t48 * mrSges(7,1) - t49 * mrSges(7,3);
t104 = -t48 * mrSges(6,1) - t49 * mrSges(6,2) - t33;
t96 = m(7) * (-t58 * pkin(5) + t47 * qJ(6) + 0.2e1 * qJD(6) * t59 - t48 * t32 + t105) + t59 * t37 + t47 * mrSges(7,3);
t11 = m(6) * t85 - t111 + (t35 + t38) * t59 + t104 * t49 + (mrSges(6,1) + mrSges(7,1)) * t47 + t106 * t30;
t51 = (mrSges(5,1) * t73 + mrSges(5,2) * t75) * qJD(1);
t57 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t100;
t9 = m(6) * t105 - t47 * mrSges(6,2) + t104 * t48 + t106 * t29 - t59 * t36 + t96;
t6 = m(5) * t101 - qJDD(4) * mrSges(5,2) + t54 * mrSges(5,3) - qJD(4) * t57 - t72 * t11 + t110 * t9 - t51 * t99;
t56 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t99;
t7 = m(5) * t90 + qJDD(4) * mrSges(5,1) - t55 * mrSges(5,3) + qJD(4) * t56 - t100 * t51 - t113;
t88 = m(4) * t69 + t75 * t6 - t73 * t7;
t86 = m(3) * t69 + t88;
t84 = -m(4) * (-qJDD(1) * pkin(2) + t82) - t73 * t6 - t75 * t7;
t81 = m(5) * t83 - t54 * mrSges(5,1) + t55 * mrSges(5,2) + t57 * t100 + t110 * t11 + t56 * t99 + t72 * t9;
t79 = -m(4) * (t78 * pkin(2) + t114) + t81;
t4 = m(3) * t102 + qJDD(1) * t108 - t109 * t78 + t79;
t3 = m(3) * t89 + qJDD(1) * t109 + t108 * t78 + t84;
t2 = m(2) * t87 - t78 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t70 * t3 + t71 * t4;
t1 = m(2) * t93 + qJDD(1) * mrSges(2,1) - t78 * mrSges(2,2) + t71 * t3 + t70 * t4;
t5 = [-m(1) * g(1) - t74 * t1 + t76 * t2, t2, t4, t88, t6, t9, -t29 * mrSges(7,2) - t48 * t33 + t96; -m(1) * g(2) + t76 * t1 + t74 * t2, t1, t3, -t78 * mrSges(4,2) - qJDD(1) * mrSges(4,3) - t79, t7, t11, -t30 * mrSges(7,3) - t49 * t37 + t95; (-m(1) - m(2)) * g(3) + t86, -m(2) * g(3) + t86, t86, qJDD(1) * mrSges(4,2) - t78 * mrSges(4,3) - t84, t81, t113, -t47 * mrSges(7,1) + t30 * mrSges(7,2) + t49 * t33 - t59 * t38 + t111;];
f_new  = t5;
