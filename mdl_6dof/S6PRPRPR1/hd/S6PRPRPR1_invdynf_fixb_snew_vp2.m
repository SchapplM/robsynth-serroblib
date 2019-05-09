% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR1
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3,theta5]';
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
% Datum: 2019-05-04 22:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:07:43
% EndTime: 2019-05-04 22:07:48
% DurationCPUTime: 2.57s
% Computational Cost: add. (35273->151), mult. (70331->211), div. (0->0), fcn. (49756->14), ass. (0->89)
t80 = sin(pkin(12));
t84 = cos(pkin(12));
t89 = sin(qJ(4));
t92 = cos(qJ(4));
t60 = (t80 * t89 - t84 * t92) * qJD(2);
t83 = sin(pkin(6));
t93 = cos(qJ(2));
t118 = t83 * t93;
t82 = sin(pkin(10));
t86 = cos(pkin(10));
t69 = g(1) * t82 - g(2) * t86;
t87 = cos(pkin(6));
t120 = t69 * t87;
t70 = -g(1) * t86 - g(2) * t82;
t79 = -g(3) + qJDD(1);
t90 = sin(qJ(2));
t104 = t118 * t79 + t120 * t93 - t70 * t90;
t41 = qJDD(2) * pkin(2) + t104;
t119 = t83 * t90;
t110 = t119 * t79 + t120 * t90 + t70 * t93;
t95 = qJD(2) ^ 2;
t42 = -pkin(2) * t95 + t110;
t81 = sin(pkin(11));
t85 = cos(pkin(11));
t106 = t85 * t41 - t42 * t81;
t100 = -qJDD(2) * pkin(3) - t106;
t113 = qJD(2) * qJD(4);
t108 = t92 * t113;
t67 = qJDD(2) * t89 + t108;
t68 = qJDD(2) * t92 - t113 * t89;
t115 = qJD(2) * t89;
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t115;
t114 = qJD(2) * t92;
t73 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t114;
t121 = 2 * qJD(5);
t116 = t41 * t81 + t42 * t85;
t32 = -pkin(3) * t95 + qJDD(2) * pkin(8) + t116;
t105 = -t69 * t83 + t79 * t87;
t55 = qJDD(3) + t105;
t107 = -t89 * t32 + t55 * t92;
t26 = (-t67 + t108) * qJ(5) + (t89 * t92 * t95 + qJDD(4)) * pkin(4) + t107;
t117 = t32 * t92 + t55 * t89;
t71 = qJD(4) * pkin(4) - qJ(5) * t115;
t78 = t92 ^ 2;
t27 = -pkin(4) * t78 * t95 + qJ(5) * t68 - qJD(4) * t71 + t117;
t112 = -t121 * t60 + t26 * t80 + t27 * t84;
t61 = (t80 * t92 + t84 * t89) * qJD(2);
t45 = pkin(5) * t60 - pkin(9) * t61;
t94 = qJD(4) ^ 2;
t22 = -pkin(5) * t94 + qJDD(4) * pkin(9) - t45 * t60 + t112;
t48 = -t67 * t80 + t68 * t84;
t49 = t67 * t84 + t68 * t80;
t97 = -t68 * pkin(4) + qJDD(5) + t71 * t115 + (-qJ(5) * t78 - pkin(8)) * t95 + t100;
t23 = (qJD(4) * t60 - t49) * pkin(9) + (qJD(4) * t61 - t48) * pkin(5) + t97;
t88 = sin(qJ(6));
t91 = cos(qJ(6));
t52 = qJD(4) * t91 - t61 * t88;
t34 = qJD(6) * t52 + qJDD(4) * t88 + t49 * t91;
t53 = qJD(4) * t88 + t61 * t91;
t35 = -mrSges(7,1) * t52 + mrSges(7,2) * t53;
t59 = qJD(6) + t60;
t36 = -mrSges(7,2) * t59 + mrSges(7,3) * t52;
t47 = qJDD(6) - t48;
t19 = m(7) * (-t22 * t88 + t23 * t91) - t34 * mrSges(7,3) + t47 * mrSges(7,1) - t53 * t35 + t59 * t36;
t33 = -qJD(6) * t53 + qJDD(4) * t91 - t49 * t88;
t37 = mrSges(7,1) * t59 - mrSges(7,3) * t53;
t20 = m(7) * (t22 * t91 + t23 * t88) + t33 * mrSges(7,3) - t47 * mrSges(7,2) + t52 * t35 - t59 * t37;
t56 = -qJD(4) * mrSges(6,2) - mrSges(6,3) * t60;
t57 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t61;
t99 = -m(6) * t97 + t48 * mrSges(6,1) - mrSges(6,2) * t49 - t19 * t91 - t20 * t88 - t60 * t56 - t57 * t61;
t122 = (t72 * t89 - t73 * t92) * qJD(2) + m(5) * (-t95 * pkin(8) + t100) - t68 * mrSges(5,1) + t67 * mrSges(5,2) - t99;
t44 = mrSges(6,1) * t60 + mrSges(6,2) * t61;
t15 = m(6) * t112 - qJDD(4) * mrSges(6,2) + mrSges(6,3) * t48 - qJD(4) * t57 - t19 * t88 + t20 * t91 - t44 * t60;
t103 = -t84 * t26 + t80 * t27;
t98 = m(7) * (-qJDD(4) * pkin(5) - t94 * pkin(9) + (t121 + t45) * t61 + t103) - t33 * mrSges(7,1) + t34 * mrSges(7,2) - t52 * t36 + t53 * t37;
t16 = m(6) * (-0.2e1 * qJD(5) * t61 - t103) - t49 * mrSges(6,3) + qJDD(4) * mrSges(6,1) - t61 * t44 + qJD(4) * t56 - t98;
t66 = (-mrSges(5,1) * t92 + mrSges(5,2) * t89) * qJD(2);
t12 = m(5) * t107 + qJDD(4) * mrSges(5,1) - t67 * mrSges(5,3) + qJD(4) * t73 - t115 * t66 + t80 * t15 + t84 * t16;
t13 = m(5) * t117 - qJDD(4) * mrSges(5,2) + mrSges(5,3) * t68 - qJD(4) * t72 + t114 * t66 + t15 * t84 - t16 * t80;
t111 = m(4) * t55 + t12 * t92 + t13 * t89;
t14 = m(4) * t106 + qJDD(2) * mrSges(4,1) - t95 * mrSges(4,2) - t122;
t7 = m(4) * t116 - mrSges(4,1) * t95 - qJDD(2) * mrSges(4,2) - t12 * t89 + t13 * t92;
t5 = m(3) * t104 + qJDD(2) * mrSges(3,1) - mrSges(3,2) * t95 + t14 * t85 + t7 * t81;
t6 = m(3) * t110 - mrSges(3,1) * t95 - qJDD(2) * mrSges(3,2) - t14 * t81 + t7 * t85;
t9 = m(3) * t105 + t111;
t109 = m(2) * t79 + t118 * t5 + t119 * t6 + t87 * t9;
t2 = m(2) * t70 - t5 * t90 + t6 * t93;
t1 = m(2) * t69 - t83 * t9 + (t5 * t93 + t6 * t90) * t87;
t3 = [-m(1) * g(1) - t1 * t82 + t2 * t86, t2, t6, t7, t13, t15, t20; -m(1) * g(2) + t1 * t86 + t2 * t82, t1, t5, t14, t12, t16, t19; -m(1) * g(3) + t109, t109, t9, t111, t122, -t99, t98;];
f_new  = t3;
