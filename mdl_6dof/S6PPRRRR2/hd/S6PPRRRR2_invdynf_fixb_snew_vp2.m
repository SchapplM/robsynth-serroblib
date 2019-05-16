% Calculate vector of cutting forces with Newton-Euler
% S6PPRRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:58
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:54:41
% EndTime: 2019-05-04 20:54:48
% DurationCPUTime: 4.71s
% Computational Cost: add. (72752->146), mult. (132608->202), div. (0->0), fcn. (103939->16), ass. (0->92)
t76 = sin(pkin(12));
t80 = cos(pkin(12));
t67 = -t80 * g(1) - t76 * g(2);
t75 = sin(pkin(13));
t79 = cos(pkin(13));
t66 = t76 * g(1) - t80 * g(2);
t74 = -g(3) + qJDD(1);
t78 = sin(pkin(6));
t82 = cos(pkin(6));
t97 = t66 * t82 + t74 * t78;
t44 = -t75 * t67 + t97 * t79;
t55 = -t78 * t66 + t82 * t74 + qJDD(2);
t77 = sin(pkin(7));
t81 = cos(pkin(7));
t117 = t44 * t81 + t55 * t77;
t89 = cos(qJ(4));
t106 = t89 * qJD(3);
t45 = t79 * t67 + t97 * t75;
t86 = sin(qJ(3));
t90 = cos(qJ(3));
t104 = t117 * t86 + t90 * t45;
t92 = qJD(3) ^ 2;
t31 = -t92 * pkin(3) + qJDD(3) * pkin(9) + t104;
t37 = -t77 * t44 + t81 * t55;
t85 = sin(qJ(4));
t108 = t89 * t31 + t85 * t37;
t63 = (-pkin(4) * t89 - pkin(10) * t85) * qJD(3);
t91 = qJD(4) ^ 2;
t24 = -t91 * pkin(4) + qJDD(4) * pkin(10) + t63 * t106 + t108;
t105 = qJD(3) * qJD(4);
t102 = t89 * t105;
t115 = t117 * t90 - t86 * t45;
t30 = -qJDD(3) * pkin(3) - t92 * pkin(9) - t115;
t64 = t85 * qJDD(3) + t102;
t72 = t85 * t105;
t65 = t89 * qJDD(3) - t72;
t27 = (-t64 - t102) * pkin(10) + (-t65 + t72) * pkin(4) + t30;
t84 = sin(qJ(5));
t88 = cos(qJ(5));
t101 = -t84 * t24 + t88 * t27;
t107 = qJD(3) * t85;
t60 = t88 * qJD(4) - t84 * t107;
t47 = t60 * qJD(5) + t84 * qJDD(4) + t88 * t64;
t59 = qJDD(5) - t65;
t61 = t84 * qJD(4) + t88 * t107;
t71 = qJD(5) - t106;
t18 = (t60 * t71 - t47) * pkin(11) + (t60 * t61 + t59) * pkin(5) + t101;
t109 = t88 * t24 + t84 * t27;
t46 = -t61 * qJD(5) + t88 * qJDD(4) - t84 * t64;
t54 = t71 * pkin(5) - t61 * pkin(11);
t58 = t60 ^ 2;
t19 = -t58 * pkin(5) + t46 * pkin(11) - t71 * t54 + t109;
t83 = sin(qJ(6));
t87 = cos(qJ(6));
t48 = t87 * t60 - t83 * t61;
t34 = t48 * qJD(6) + t83 * t46 + t87 * t47;
t49 = t83 * t60 + t87 * t61;
t38 = -t48 * mrSges(7,1) + t49 * mrSges(7,2);
t70 = qJD(6) + t71;
t40 = -t70 * mrSges(7,2) + t48 * mrSges(7,3);
t57 = qJDD(6) + t59;
t16 = m(7) * (t87 * t18 - t83 * t19) - t34 * mrSges(7,3) + t57 * mrSges(7,1) - t49 * t38 + t70 * t40;
t33 = -t49 * qJD(6) + t87 * t46 - t83 * t47;
t41 = t70 * mrSges(7,1) - t49 * mrSges(7,3);
t17 = m(7) * (t83 * t18 + t87 * t19) + t33 * mrSges(7,3) - t57 * mrSges(7,2) + t48 * t38 - t70 * t41;
t50 = -t60 * mrSges(6,1) + t61 * mrSges(6,2);
t52 = -t71 * mrSges(6,2) + t60 * mrSges(6,3);
t13 = m(6) * t101 + t59 * mrSges(6,1) - t47 * mrSges(6,3) + t87 * t16 + t83 * t17 - t61 * t50 + t71 * t52;
t53 = t71 * mrSges(6,1) - t61 * mrSges(6,3);
t14 = m(6) * t109 - t59 * mrSges(6,2) + t46 * mrSges(6,3) - t83 * t16 + t87 * t17 + t60 * t50 - t71 * t53;
t62 = (-mrSges(5,1) * t89 + mrSges(5,2) * t85) * qJD(3);
t68 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t107;
t12 = m(5) * t108 - qJDD(4) * mrSges(5,2) + t65 * mrSges(5,3) - qJD(4) * t68 + t62 * t106 - t84 * t13 + t88 * t14;
t100 = -t85 * t31 + t89 * t37;
t69 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t106;
t23 = -qJDD(4) * pkin(4) - t91 * pkin(10) + t63 * t107 - t100;
t95 = t33 * mrSges(7,1) + t48 * t40 - m(7) * (-t46 * pkin(5) - t58 * pkin(11) + t61 * t54 + t23) - t34 * mrSges(7,2) - t49 * t41;
t93 = m(6) * t23 - t46 * mrSges(6,1) + t47 * mrSges(6,2) - t60 * t52 + t61 * t53 - t95;
t15 = m(5) * t100 + qJDD(4) * mrSges(5,1) - t64 * mrSges(5,3) + qJD(4) * t69 - t62 * t107 - t93;
t10 = m(4) * t37 + t85 * t12 + t89 * t15;
t114 = m(5) * t30 - t65 * mrSges(5,1) + t64 * mrSges(5,2) + t88 * t13 + t84 * t14 + (t85 * t68 - t89 * t69) * qJD(3);
t11 = m(4) * t115 + qJDD(3) * mrSges(4,1) - t92 * mrSges(4,2) - t114;
t9 = m(4) * t104 - t92 * mrSges(4,1) - qJDD(3) * mrSges(4,2) + t89 * t12 - t85 * t15;
t99 = t11 * t90 + t86 * t9;
t4 = m(3) * t44 - t77 * t10 + t99 * t81;
t8 = m(3) * t45 - t86 * t11 + t90 * t9;
t116 = t4 * t79 + t75 * t8;
t6 = m(3) * t55 + t81 * t10 + t99 * t77;
t103 = m(2) * t74 + t116 * t78 + t82 * t6;
t2 = m(2) * t67 - t75 * t4 + t79 * t8;
t1 = m(2) * t66 + t116 * t82 - t78 * t6;
t3 = [-m(1) * g(1) - t76 * t1 + t80 * t2, t2, t8, t9, t12, t14, t17; -m(1) * g(2) + t80 * t1 + t76 * t2, t1, t4, t11, t15, t13, t16; -m(1) * g(3) + t103, t103, t6, t10, t114, t93, -t95;];
f_new  = t3;
