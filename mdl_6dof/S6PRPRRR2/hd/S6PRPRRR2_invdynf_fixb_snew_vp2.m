% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1,theta3]';
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
% Datum: 2019-05-05 00:26
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 00:21:42
% EndTime: 2019-05-05 00:21:48
% DurationCPUTime: 2.94s
% Computational Cost: add. (43418->150), mult. (80752->204), div. (0->0), fcn. (57411->14), ass. (0->91)
t93 = cos(qJ(4));
t111 = t93 * qJD(2);
t83 = sin(pkin(6));
t94 = cos(qJ(2));
t116 = t83 * t94;
t82 = sin(pkin(11));
t85 = cos(pkin(11));
t69 = t82 * g(1) - t85 * g(2);
t86 = cos(pkin(6));
t118 = t69 * t86;
t70 = -t85 * g(1) - t82 * g(2);
t80 = -g(3) + qJDD(1);
t90 = sin(qJ(2));
t101 = t80 * t116 + t94 * t118 - t90 * t70;
t43 = qJDD(2) * pkin(2) + t101;
t117 = t83 * t90;
t108 = t80 * t117 + t90 * t118 + t94 * t70;
t96 = qJD(2) ^ 2;
t44 = -t96 * pkin(2) + t108;
t81 = sin(pkin(12));
t84 = cos(pkin(12));
t113 = t81 * t43 + t84 * t44;
t35 = -t96 * pkin(3) + qJDD(2) * pkin(8) + t113;
t102 = -t83 * t69 + t86 * t80;
t55 = qJDD(3) + t102;
t89 = sin(qJ(4));
t114 = t93 * t35 + t89 * t55;
t66 = (-pkin(4) * t93 - pkin(9) * t89) * qJD(2);
t95 = qJD(4) ^ 2;
t25 = -t95 * pkin(4) + qJDD(4) * pkin(9) + t66 * t111 + t114;
t110 = qJD(2) * qJD(4);
t106 = t93 * t110;
t103 = t84 * t43 - t81 * t44;
t34 = -qJDD(2) * pkin(3) - t96 * pkin(8) - t103;
t67 = t89 * qJDD(2) + t106;
t78 = t89 * t110;
t68 = t93 * qJDD(2) - t78;
t28 = (-t67 - t106) * pkin(9) + (-t68 + t78) * pkin(4) + t34;
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t105 = -t88 * t25 + t92 * t28;
t112 = qJD(2) * t89;
t63 = t92 * qJD(4) - t88 * t112;
t46 = t63 * qJD(5) + t88 * qJDD(4) + t92 * t67;
t61 = qJDD(5) - t68;
t64 = t88 * qJD(4) + t92 * t112;
t76 = qJD(5) - t111;
t19 = (t63 * t76 - t46) * pkin(10) + (t63 * t64 + t61) * pkin(5) + t105;
t115 = t92 * t25 + t88 * t28;
t45 = -t64 * qJD(5) + t92 * qJDD(4) - t88 * t67;
t54 = t76 * pkin(5) - t64 * pkin(10);
t60 = t63 ^ 2;
t20 = -t60 * pkin(5) + t45 * pkin(10) - t76 * t54 + t115;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t47 = t91 * t63 - t87 * t64;
t31 = t47 * qJD(6) + t87 * t45 + t91 * t46;
t48 = t87 * t63 + t91 * t64;
t37 = -t47 * mrSges(7,1) + t48 * mrSges(7,2);
t75 = qJD(6) + t76;
t41 = -t75 * mrSges(7,2) + t47 * mrSges(7,3);
t59 = qJDD(6) + t61;
t17 = m(7) * (t91 * t19 - t87 * t20) - t31 * mrSges(7,3) + t59 * mrSges(7,1) - t48 * t37 + t75 * t41;
t30 = -t48 * qJD(6) + t91 * t45 - t87 * t46;
t42 = t75 * mrSges(7,1) - t48 * mrSges(7,3);
t18 = m(7) * (t87 * t19 + t91 * t20) + t30 * mrSges(7,3) - t59 * mrSges(7,2) + t47 * t37 - t75 * t42;
t49 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t52 = -t76 * mrSges(6,2) + t63 * mrSges(6,3);
t13 = m(6) * t105 + t61 * mrSges(6,1) - t46 * mrSges(6,3) + t91 * t17 + t87 * t18 - t64 * t49 + t76 * t52;
t53 = t76 * mrSges(6,1) - t64 * mrSges(6,3);
t14 = m(6) * t115 - t61 * mrSges(6,2) + t45 * mrSges(6,3) - t87 * t17 + t91 * t18 + t63 * t49 - t76 * t53;
t71 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t112;
t72 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t111;
t119 = m(5) * t34 - t68 * mrSges(5,1) + t67 * mrSges(5,2) + t92 * t13 + t88 * t14 + (t89 * t71 - t93 * t72) * qJD(2);
t65 = (-mrSges(5,1) * t93 + mrSges(5,2) * t89) * qJD(2);
t12 = m(5) * t114 - qJDD(4) * mrSges(5,2) + t68 * mrSges(5,3) - qJD(4) * t71 + t65 * t111 - t88 * t13 + t92 * t14;
t104 = -t89 * t35 + t93 * t55;
t24 = -qJDD(4) * pkin(4) - t95 * pkin(9) + t66 * t112 - t104;
t99 = t30 * mrSges(7,1) + t47 * t41 - m(7) * (-t45 * pkin(5) - t60 * pkin(10) + t64 * t54 + t24) - t31 * mrSges(7,2) - t48 * t42;
t97 = m(6) * t24 - t45 * mrSges(6,1) + t46 * mrSges(6,2) - t63 * t52 + t64 * t53 - t99;
t16 = m(5) * t104 + qJDD(4) * mrSges(5,1) - t67 * mrSges(5,3) + qJD(4) * t72 - t65 * t112 - t97;
t109 = m(4) * t55 + t89 * t12 + t93 * t16;
t10 = m(4) * t103 + qJDD(2) * mrSges(4,1) - t96 * mrSges(4,2) - t119;
t7 = m(4) * t113 - t96 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t93 * t12 - t89 * t16;
t5 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t96 * mrSges(3,2) + t84 * t10 + t81 * t7;
t6 = m(3) * t108 - t96 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t81 * t10 + t84 * t7;
t9 = m(3) * t102 + t109;
t107 = m(2) * t80 + t5 * t116 + t6 * t117 + t86 * t9;
t2 = m(2) * t70 - t90 * t5 + t94 * t6;
t1 = m(2) * t69 - t83 * t9 + (t5 * t94 + t6 * t90) * t86;
t3 = [-m(1) * g(1) - t82 * t1 + t85 * t2, t2, t6, t7, t12, t14, t18; -m(1) * g(2) + t85 * t1 + t82 * t2, t1, t5, t10, t16, t13, t17; -m(1) * g(3) + t107, t107, t9, t109, t119, t97, -t99;];
f_new  = t3;
