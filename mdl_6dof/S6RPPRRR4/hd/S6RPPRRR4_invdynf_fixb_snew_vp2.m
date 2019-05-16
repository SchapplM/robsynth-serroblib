% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 15:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:40:12
% EndTime: 2019-05-05 15:40:16
% DurationCPUTime: 1.43s
% Computational Cost: add. (20659->154), mult. (37570->194), div. (0->0), fcn. (20499->10), ass. (0->83)
t81 = qJD(1) ^ 2;
t104 = -pkin(1) - pkin(2);
t75 = sin(qJ(1));
t79 = cos(qJ(1));
t91 = -t79 * g(1) - t75 * g(2);
t88 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t91;
t44 = t104 * t81 + t88;
t97 = t75 * g(1) - t79 * g(2);
t85 = -t81 * qJ(2) + qJDD(2) - t97;
t46 = t104 * qJDD(1) + t85;
t70 = sin(pkin(10));
t71 = cos(pkin(10));
t92 = -t70 * t44 + t71 * t46;
t29 = qJDD(1) * pkin(3) - t81 * pkin(7) - t92;
t74 = sin(qJ(4));
t78 = cos(qJ(4));
t98 = qJD(1) * qJD(4);
t95 = t78 * t98;
t55 = -t74 * qJDD(1) - t95;
t96 = t74 * t98;
t56 = -t78 * qJDD(1) + t96;
t99 = qJD(1) * t74;
t57 = (qJD(4) * mrSges(5,1)) + mrSges(5,3) * t99;
t63 = t78 * qJD(1);
t58 = -(qJD(4) * mrSges(5,2)) - mrSges(5,3) * t63;
t73 = sin(qJ(5));
t77 = cos(qJ(5));
t51 = t77 * qJD(4) + t73 * t99;
t34 = t51 * qJD(5) + t73 * qJDD(4) + t77 * t55;
t50 = qJDD(5) - t56;
t52 = t73 * qJD(4) - t77 * t99;
t60 = t63 + qJD(5);
t21 = (-t55 + t95) * pkin(8) + (-t56 - t96) * pkin(4) + t29;
t100 = t71 * t44 + t70 * t46;
t30 = -t81 * pkin(3) - qJDD(1) * pkin(7) + t100;
t68 = g(3) + qJDD(3);
t101 = t78 * t30 + t74 * t68;
t54 = (pkin(4) * t78 + pkin(8) * t74) * qJD(1);
t80 = qJD(4) ^ 2;
t24 = -t80 * pkin(4) + qJDD(4) * pkin(8) - t54 * t63 + t101;
t94 = t77 * t21 - t73 * t24;
t12 = (t51 * t60 - t34) * pkin(9) + (t51 * t52 + t50) * pkin(5) + t94;
t102 = t73 * t21 + t77 * t24;
t33 = -t52 * qJD(5) + t77 * qJDD(4) - t73 * t55;
t43 = t60 * pkin(5) - t52 * pkin(9);
t49 = t51 ^ 2;
t13 = -t49 * pkin(5) + t33 * pkin(9) - t60 * t43 + t102;
t72 = sin(qJ(6));
t76 = cos(qJ(6));
t35 = t76 * t51 - t72 * t52;
t18 = t35 * qJD(6) + t72 * t33 + t76 * t34;
t36 = t72 * t51 + t76 * t52;
t26 = -t35 * mrSges(7,1) + t36 * mrSges(7,2);
t59 = qJD(6) + t60;
t31 = -t59 * mrSges(7,2) + t35 * mrSges(7,3);
t48 = qJDD(6) + t50;
t10 = m(7) * (t76 * t12 - t72 * t13) - t18 * mrSges(7,3) + t48 * mrSges(7,1) - t36 * t26 + t59 * t31;
t17 = -t36 * qJD(6) + t76 * t33 - t72 * t34;
t32 = t59 * mrSges(7,1) - t36 * mrSges(7,3);
t11 = m(7) * (t72 * t12 + t76 * t13) + t17 * mrSges(7,3) - t48 * mrSges(7,2) + t35 * t26 - t59 * t32;
t37 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t41 = -t60 * mrSges(6,2) + t51 * mrSges(6,3);
t7 = m(6) * t94 + t50 * mrSges(6,1) - t34 * mrSges(6,3) + t76 * t10 + t72 * t11 - t52 * t37 + t60 * t41;
t42 = t60 * mrSges(6,1) - t52 * mrSges(6,3);
t8 = m(6) * t102 - t50 * mrSges(6,2) + t33 * mrSges(6,3) - t72 * t10 + t76 * t11 + t51 * t37 - t60 * t42;
t106 = m(5) * t29 - t56 * mrSges(5,1) + t55 * mrSges(5,2) - (t74 * t57 - t78 * t58) * qJD(1) + t77 * t7 + t73 * t8;
t105 = -m(2) - m(3);
t103 = mrSges(2,1) + mrSges(3,1);
t93 = -t74 * t30 + t78 * t68;
t53 = (mrSges(5,1) * t78 - mrSges(5,2) * t74) * qJD(1);
t6 = m(5) * t101 - qJDD(4) * mrSges(5,2) + t56 * mrSges(5,3) - qJD(4) * t57 - t53 * t63 - t73 * t7 + t77 * t8;
t23 = -qJDD(4) * pkin(4) - t80 * pkin(8) - t54 * t99 - t93;
t84 = t17 * mrSges(7,1) + t35 * t31 - m(7) * (-t33 * pkin(5) - t49 * pkin(9) + t52 * t43 + t23) - t18 * mrSges(7,2) - t36 * t32;
t82 = m(6) * t23 - t33 * mrSges(6,1) + t34 * mrSges(6,2) - t51 * t41 + t52 * t42 - t84;
t9 = m(5) * t93 + qJDD(4) * mrSges(5,1) - t55 * mrSges(5,3) + qJD(4) * t58 + t53 * t99 - t82;
t4 = m(4) * t100 - t81 * mrSges(4,1) + qJDD(1) * mrSges(4,2) + t78 * t6 - t74 * t9;
t5 = m(4) * t92 - qJDD(1) * mrSges(4,1) - t81 * mrSges(4,2) - t106;
t89 = t71 * t4 - t70 * t5 + m(3) * (-t81 * pkin(1) + t88) + qJDD(1) * mrSges(3,3);
t87 = -m(3) * (-qJDD(1) * pkin(1) + t85) - t70 * t4 - t71 * t5;
t86 = m(4) * t68 + t74 * t6 + t78 * t9;
t2 = m(2) * t97 + (-mrSges(2,2) + mrSges(3,3)) * t81 + t103 * qJDD(1) + t87;
t1 = m(2) * t91 - qJDD(1) * mrSges(2,2) - t103 * t81 + t89;
t3 = [-m(1) * g(1) + t79 * t1 - t75 * t2, t1, -t81 * mrSges(3,1) + t89, t4, t6, t8, t11; -m(1) * g(2) + t75 * t1 + t79 * t2, t2, -m(3) * g(3) - t86, t5, t9, t7, t10; (-m(1) + t105) * g(3) - t86, t105 * g(3) - t86, -qJDD(1) * mrSges(3,1) - t81 * mrSges(3,3) - t87, t86, t106, t82, -t84;];
f_new  = t3;
