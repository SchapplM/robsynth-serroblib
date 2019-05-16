% Calculate vector of cutting forces with Newton-Euler
% S6RPPPRR3
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
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-05-05 13:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPPRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:40:45
% EndTime: 2019-05-05 13:40:47
% DurationCPUTime: 1.18s
% Computational Cost: add. (15247->142), mult. (30775->178), div. (0->0), fcn. (17601->10), ass. (0->81)
t77 = qJD(1) ^ 2;
t68 = cos(pkin(10));
t61 = t68 ^ 2;
t66 = sin(pkin(10));
t101 = t66 ^ 2 + t61;
t109 = t101 * mrSges(5,3);
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t89 = t66 * t71 - t68 * t74;
t50 = t89 * qJD(1);
t90 = -t66 * t74 - t68 * t71;
t51 = t90 * qJD(1);
t99 = t51 * qJD(5);
t36 = t89 * qJDD(1) - t99;
t108 = -m(2) - m(3);
t107 = -pkin(1) - pkin(2);
t106 = pkin(4) * t77;
t105 = mrSges(2,1) + mrSges(3,1);
t64 = g(3) + qJDD(3);
t97 = qJD(1) * qJD(4);
t102 = t68 * t64 + 0.2e1 * t66 * t97;
t72 = sin(qJ(1));
t75 = cos(qJ(1));
t93 = -t75 * g(1) - t72 * g(2);
t85 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t93;
t43 = t107 * t77 + t85;
t95 = t72 * g(1) - t75 * g(2);
t82 = -t77 * qJ(2) + qJDD(2) - t95;
t47 = t107 * qJDD(1) + t82;
t67 = sin(pkin(9));
t69 = cos(pkin(9));
t103 = t69 * t43 + t67 * t47;
t27 = -t77 * pkin(3) - qJDD(1) * qJ(4) + t103;
t19 = (pkin(7) * qJDD(1) + t68 * t106 - t27) * t66 + t102;
t96 = t66 * t64 + (t27 - 0.2e1 * t97) * t68;
t98 = qJDD(1) * t68;
t20 = -pkin(7) * t98 - t61 * t106 + t96;
t104 = t71 * t19 + t74 * t20;
t100 = t50 * qJD(5);
t94 = -t67 * t43 + t69 * t47;
t92 = t68 * mrSges(5,1) - t66 * mrSges(5,2);
t91 = t74 * t19 - t71 * t20;
t35 = -t50 * pkin(5) - t51 * pkin(8);
t76 = qJD(5) ^ 2;
t15 = -t76 * pkin(5) + qJDD(5) * pkin(8) + t50 * t35 + t104;
t37 = t90 * qJDD(1) + t100;
t86 = qJDD(1) * pkin(3) + qJDD(4) - t94;
t78 = pkin(4) * t98 + (-t101 * pkin(7) - qJ(4)) * t77 + t86;
t16 = (-t37 - t100) * pkin(8) + (-t36 + t99) * pkin(5) + t78;
t70 = sin(qJ(6));
t73 = cos(qJ(6));
t40 = t73 * qJD(5) - t70 * t51;
t24 = t40 * qJD(6) + t70 * qJDD(5) + t73 * t37;
t41 = t70 * qJD(5) + t73 * t51;
t28 = -t40 * mrSges(7,1) + t41 * mrSges(7,2);
t48 = qJD(6) - t50;
t29 = -t48 * mrSges(7,2) + t40 * mrSges(7,3);
t34 = qJDD(6) - t36;
t12 = m(7) * (-t70 * t15 + t73 * t16) - t24 * mrSges(7,3) + t34 * mrSges(7,1) - t41 * t28 + t48 * t29;
t23 = -t41 * qJD(6) + t73 * qJDD(5) - t70 * t37;
t30 = t48 * mrSges(7,1) - t41 * mrSges(7,3);
t13 = m(7) * (t73 * t15 + t70 * t16) + t23 * mrSges(7,3) - t34 * mrSges(7,2) + t40 * t28 - t48 * t30;
t32 = -t50 * mrSges(6,1) + t51 * mrSges(6,2);
t45 = qJD(5) * mrSges(6,1) - t51 * mrSges(6,3);
t8 = m(6) * t104 - qJDD(5) * mrSges(6,2) + t36 * mrSges(6,3) - qJD(5) * t45 - t70 * t12 + t73 * t13 + t50 * t32;
t87 = -qJDD(1) * mrSges(5,3) - t77 * t92;
t44 = -qJD(5) * mrSges(6,2) + t50 * mrSges(6,3);
t79 = m(7) * (-qJDD(5) * pkin(5) - t76 * pkin(8) + t51 * t35 - t91) - t23 * mrSges(7,1) + t24 * mrSges(7,2) - t40 * t29 + t41 * t30;
t9 = m(6) * t91 + qJDD(5) * mrSges(6,1) - t37 * mrSges(6,3) + qJD(5) * t44 - t51 * t32 - t79;
t5 = m(5) * t102 + t71 * t8 + t74 * t9 + (-m(5) * t27 - t87) * t66;
t6 = m(5) * t96 + t87 * t68 - t71 * t9 + t74 * t8;
t4 = m(4) * t103 - t77 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t66 * t5 + t68 * t6;
t81 = -m(6) * t78 + t36 * mrSges(6,1) - t37 * mrSges(6,2) - t73 * t12 - t70 * t13 + t50 * t44 - t51 * t45;
t80 = m(5) * (-t77 * qJ(4) + t86) - t81;
t7 = m(4) * t94 + (-mrSges(4,2) + t109) * t77 + (-mrSges(4,1) - t92) * qJDD(1) - t80;
t88 = t69 * t4 - t67 * t7 + m(3) * (-t77 * pkin(1) + t85) + qJDD(1) * mrSges(3,3);
t84 = -m(3) * (-qJDD(1) * pkin(1) + t82) - t67 * t4 - t69 * t7;
t83 = m(4) * t64 + t68 * t5 + t66 * t6;
t2 = m(2) * t95 + (-mrSges(2,2) + mrSges(3,3)) * t77 + t105 * qJDD(1) + t84;
t1 = m(2) * t93 - qJDD(1) * mrSges(2,2) - t105 * t77 + t88;
t3 = [-m(1) * g(1) + t75 * t1 - t72 * t2, t1, -t77 * mrSges(3,1) + t88, t4, t6, t8, t13; -m(1) * g(2) + t72 * t1 + t75 * t2, t2, -m(3) * g(3) - t83, t7, t5, t9, t12; (-m(1) + t108) * g(3) - t83, t108 * g(3) - t83, -qJDD(1) * mrSges(3,1) - t77 * mrSges(3,3) - t84, t83, t92 * qJDD(1) - t77 * t109 + t80, -t81, t79;];
f_new  = t3;
