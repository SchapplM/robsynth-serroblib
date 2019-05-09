% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-05-05 15:32
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:30:50
% EndTime: 2019-05-05 15:30:53
% DurationCPUTime: 1.30s
% Computational Cost: add. (16942->151), mult. (31705->191), div. (0->0), fcn. (18986->10), ass. (0->84)
t83 = sin(qJ(1));
t87 = cos(qJ(1));
t104 = t83 * g(1) - t87 * g(2);
t57 = qJDD(1) * pkin(1) + t104;
t89 = qJD(1) ^ 2;
t98 = -t87 * g(1) - t83 * g(2);
t59 = -t89 * pkin(1) + t98;
t78 = sin(pkin(10));
t79 = cos(pkin(10));
t110 = t78 * t57 + t79 * t59;
t115 = -qJDD(1) * qJ(3) - (2 * qJD(3) * qJD(1)) - t110;
t114 = -pkin(2) - pkin(7);
t113 = mrSges(3,1) - mrSges(4,2);
t112 = -mrSges(3,2) + mrSges(4,3);
t107 = qJD(1) * qJD(4);
t82 = sin(qJ(4));
t103 = t82 * t107;
t86 = cos(qJ(4));
t69 = t86 * t107;
t61 = -t82 * qJDD(1) - t69;
t62 = t86 * qJDD(1) - t103;
t95 = t114 * t89 - t115;
t23 = (-t62 + t103) * pkin(8) + (-t61 + t69) * pkin(4) + t95;
t100 = t79 * t57 - t78 * t59;
t94 = -t89 * qJ(3) + qJDD(3) - t100;
t34 = t114 * qJDD(1) + t94;
t77 = -g(3) + qJDD(2);
t109 = t82 * t34 + t86 * t77;
t60 = (pkin(4) * t82 - pkin(8) * t86) * qJD(1);
t71 = t82 * qJD(1);
t88 = qJD(4) ^ 2;
t27 = -t88 * pkin(4) + qJDD(4) * pkin(8) - t60 * t71 + t109;
t81 = sin(qJ(5));
t85 = cos(qJ(5));
t111 = t81 * t23 + t85 * t27;
t108 = qJD(1) * t86;
t66 = t71 + qJD(5);
t54 = qJDD(5) - t61;
t102 = t85 * t23 - t81 * t27;
t101 = t86 * t34 - t82 * t77;
t58 = (mrSges(5,1) * t82 + mrSges(5,2) * t86) * qJD(1);
t63 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t71;
t26 = -qJDD(4) * pkin(4) - t88 * pkin(8) + t60 * t108 - t101;
t56 = t81 * qJD(4) + t85 * t108;
t37 = -t56 * qJD(5) + t85 * qJDD(4) - t81 * t62;
t55 = t85 * qJD(4) - t81 * t108;
t38 = t55 * qJD(5) + t81 * qJDD(4) + t85 * t62;
t42 = -t66 * mrSges(6,2) + t55 * mrSges(6,3);
t43 = t66 * mrSges(6,1) - t56 * mrSges(6,3);
t80 = sin(qJ(6));
t84 = cos(qJ(6));
t40 = t80 * t55 + t84 * t56;
t19 = -t40 * qJD(6) + t84 * t37 - t80 * t38;
t39 = t84 * t55 - t80 * t56;
t20 = t39 * qJD(6) + t80 * t37 + t84 * t38;
t65 = qJD(6) + t66;
t32 = -t65 * mrSges(7,2) + t39 * mrSges(7,3);
t33 = t65 * mrSges(7,1) - t40 * mrSges(7,3);
t44 = t66 * pkin(5) - t56 * pkin(9);
t53 = t55 ^ 2;
t93 = t19 * mrSges(7,1) + t39 * t32 - m(7) * (-t37 * pkin(5) - t53 * pkin(9) + t56 * t44 + t26) - t20 * mrSges(7,2) - t40 * t33;
t90 = m(6) * t26 - t37 * mrSges(6,1) + t38 * mrSges(6,2) - t55 * t42 + t56 * t43 - t93;
t11 = m(5) * t101 + qJDD(4) * mrSges(5,1) - t62 * mrSges(5,3) + qJD(4) * t63 - t58 * t108 - t90;
t14 = (t55 * t66 - t38) * pkin(9) + (t55 * t56 + t54) * pkin(5) + t102;
t15 = -t53 * pkin(5) + t37 * pkin(9) - t66 * t44 + t111;
t28 = -t39 * mrSges(7,1) + t40 * mrSges(7,2);
t49 = qJDD(6) + t54;
t12 = m(7) * (t84 * t14 - t80 * t15) - t20 * mrSges(7,3) + t49 * mrSges(7,1) - t40 * t28 + t65 * t32;
t13 = m(7) * (t80 * t14 + t84 * t15) + t19 * mrSges(7,3) - t49 * mrSges(7,2) + t39 * t28 - t65 * t33;
t41 = -t55 * mrSges(6,1) + t56 * mrSges(6,2);
t10 = m(6) * t111 - t54 * mrSges(6,2) + t37 * mrSges(6,3) - t80 * t12 + t84 * t13 + t55 * t41 - t66 * t43;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108;
t9 = m(6) * t102 + t54 * mrSges(6,1) - t38 * mrSges(6,3) + t84 * t12 + t80 * t13 - t56 * t41 + t66 * t42;
t6 = m(5) * t109 - qJDD(4) * mrSges(5,2) + t61 * mrSges(5,3) - qJD(4) * t64 + t85 * t10 - t58 * t71 - t81 * t9;
t99 = m(4) * t77 - t82 * t11 + t86 * t6;
t97 = m(3) * t77 + t99;
t96 = -m(4) * (-qJDD(1) * pkin(2) + t94) - t86 * t11 - t82 * t6;
t92 = m(5) * t95 - t61 * mrSges(5,1) + t62 * mrSges(5,2) + t81 * t10 + t64 * t108 + t63 * t71 + t85 * t9;
t91 = -m(4) * (t89 * pkin(2) + t115) + t92;
t4 = m(3) * t110 + t112 * qJDD(1) - t113 * t89 + t91;
t3 = m(3) * t100 + t113 * qJDD(1) + t112 * t89 + t96;
t2 = m(2) * t98 - t89 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t78 * t3 + t79 * t4;
t1 = m(2) * t104 + qJDD(1) * mrSges(2,1) - t89 * mrSges(2,2) + t79 * t3 + t78 * t4;
t5 = [-m(1) * g(1) - t83 * t1 + t87 * t2, t2, t4, t99, t6, t10, t13; -m(1) * g(2) + t87 * t1 + t83 * t2, t1, t3, -t89 * mrSges(4,2) - qJDD(1) * mrSges(4,3) - t91, t11, t9, t12; (-m(1) - m(2)) * g(3) + t97, -m(2) * g(3) + t97, t97, qJDD(1) * mrSges(4,2) - t89 * mrSges(4,3) - t96, t92, t90, -t93;];
f_new  = t5;
