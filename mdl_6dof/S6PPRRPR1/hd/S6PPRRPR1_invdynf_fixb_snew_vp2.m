% Calculate vector of cutting forces with Newton-Euler
% S6PPRRPR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d6,theta1,theta2,theta5]';
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
% Datum: 2019-05-04 20:07
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRPR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRPR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:03:20
% EndTime: 2019-05-04 20:03:25
% DurationCPUTime: 4.60s
% Computational Cost: add. (68537->145), mult. (129854->204), div. (0->0), fcn. (101353->16), ass. (0->90)
t77 = sin(pkin(11));
t82 = cos(pkin(11));
t68 = -t82 * g(1) - t77 * g(2);
t76 = sin(pkin(12));
t81 = cos(pkin(12));
t67 = t77 * g(1) - t82 * g(2);
t74 = -g(3) + qJDD(1);
t79 = sin(pkin(6));
t84 = cos(pkin(6));
t97 = t67 * t84 + t74 * t79;
t44 = -t76 * t68 + t97 * t81;
t50 = -t79 * t67 + t84 * t74 + qJDD(2);
t78 = sin(pkin(7));
t83 = cos(pkin(7));
t117 = t44 * t83 + t50 * t78;
t89 = cos(qJ(4));
t107 = t89 * qJD(3);
t45 = t81 * t68 + t97 * t76;
t87 = sin(qJ(3));
t90 = cos(qJ(3));
t104 = t117 * t87 + t90 * t45;
t92 = qJD(3) ^ 2;
t31 = -t92 * pkin(3) + qJDD(3) * pkin(9) + t104;
t38 = -t78 * t44 + t83 * t50;
t86 = sin(qJ(4));
t109 = t89 * t31 + t86 * t38;
t63 = (-pkin(4) * t89 - qJ(5) * t86) * qJD(3);
t91 = qJD(4) ^ 2;
t24 = -t91 * pkin(4) + qJDD(4) * qJ(5) + t63 * t107 + t109;
t106 = qJD(3) * qJD(4);
t102 = t89 * t106;
t115 = t117 * t90 - t87 * t45;
t30 = -qJDD(3) * pkin(3) - t92 * pkin(9) - t115;
t65 = t86 * qJDD(3) + t102;
t72 = t86 * t106;
t66 = t89 * qJDD(3) - t72;
t27 = (-t65 - t102) * qJ(5) + (-t66 + t72) * pkin(4) + t30;
t108 = qJD(3) * t86;
t75 = sin(pkin(13));
t80 = cos(pkin(13));
t61 = t75 * qJD(4) + t80 * t108;
t100 = -0.2e1 * qJD(5) * t61 - t75 * t24 + t80 * t27;
t54 = t75 * qJDD(4) + t80 * t65;
t60 = t80 * qJD(4) - t75 * t108;
t18 = (-t60 * t107 - t54) * pkin(10) + (t60 * t61 - t66) * pkin(5) + t100;
t105 = 0.2e1 * qJD(5) * t60 + t80 * t24 + t75 * t27;
t53 = t80 * qJDD(4) - t75 * t65;
t55 = -pkin(5) * t107 - t61 * pkin(10);
t59 = t60 ^ 2;
t19 = -t59 * pkin(5) + t53 * pkin(10) + t55 * t107 + t105;
t85 = sin(qJ(6));
t88 = cos(qJ(6));
t46 = t88 * t60 - t85 * t61;
t34 = t46 * qJD(6) + t85 * t53 + t88 * t54;
t47 = t85 * t60 + t88 * t61;
t37 = -t46 * mrSges(7,1) + t47 * mrSges(7,2);
t71 = qJD(6) - t107;
t40 = -t71 * mrSges(7,2) + t46 * mrSges(7,3);
t62 = qJDD(6) - t66;
t15 = m(7) * (t88 * t18 - t85 * t19) - t34 * mrSges(7,3) + t62 * mrSges(7,1) - t47 * t37 + t71 * t40;
t33 = -t47 * qJD(6) + t88 * t53 - t85 * t54;
t41 = t71 * mrSges(7,1) - t47 * mrSges(7,3);
t16 = m(7) * (t85 * t18 + t88 * t19) + t33 * mrSges(7,3) - t62 * mrSges(7,2) + t46 * t37 - t71 * t41;
t48 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t51 = mrSges(6,2) * t107 + t60 * mrSges(6,3);
t13 = m(6) * t100 - t66 * mrSges(6,1) - t54 * mrSges(6,3) - t51 * t107 + t88 * t15 + t85 * t16 - t61 * t48;
t52 = -mrSges(6,1) * t107 - t61 * mrSges(6,3);
t14 = m(6) * t105 + t66 * mrSges(6,2) + t53 * mrSges(6,3) + t52 * t107 - t85 * t15 + t88 * t16 + t60 * t48;
t64 = (-mrSges(5,1) * t89 + mrSges(5,2) * t86) * qJD(3);
t69 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t108;
t12 = m(5) * t109 - qJDD(4) * mrSges(5,2) + t66 * mrSges(5,3) - qJD(4) * t69 + t64 * t107 - t75 * t13 + t80 * t14;
t101 = -t86 * t31 + t89 * t38;
t70 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t107;
t23 = -qJDD(4) * pkin(4) - t91 * qJ(5) + t63 * t108 + qJDD(5) - t101;
t95 = t33 * mrSges(7,1) + t46 * t40 - m(7) * (-t53 * pkin(5) - t59 * pkin(10) + t61 * t55 + t23) - t34 * mrSges(7,2) - t47 * t41;
t93 = m(6) * t23 - t53 * mrSges(6,1) + t54 * mrSges(6,2) - t60 * t51 + t61 * t52 - t95;
t17 = m(5) * t101 + qJDD(4) * mrSges(5,1) - t65 * mrSges(5,3) + qJD(4) * t70 - t64 * t108 - t93;
t10 = m(4) * t38 + t86 * t12 + t89 * t17;
t114 = m(5) * t30 - t66 * mrSges(5,1) + t65 * mrSges(5,2) + t80 * t13 + t75 * t14 + (t86 * t69 - t89 * t70) * qJD(3);
t11 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t92 * mrSges(4,2) - t114;
t9 = m(4) * t104 - t92 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t89 * t12 - t86 * t17;
t99 = t11 * t90 + t87 * t9;
t4 = m(3) * t44 - t78 * t10 + t99 * t83;
t8 = m(3) * t45 - t87 * t11 + t90 * t9;
t116 = t4 * t81 + t76 * t8;
t6 = m(3) * t50 + t83 * t10 + t99 * t78;
t103 = m(2) * t74 + t116 * t79 + t84 * t6;
t2 = m(2) * t68 - t76 * t4 + t81 * t8;
t1 = m(2) * t67 + t116 * t84 - t79 * t6;
t3 = [-m(1) * g(1) - t77 * t1 + t82 * t2, t2, t8, t9, t12, t14, t16; -m(1) * g(2) + t82 * t1 + t77 * t2, t1, t4, t11, t17, t13, t15; -m(1) * g(3) + t103, t103, t6, t10, t114, t93, -t95;];
f_new  = t3;
