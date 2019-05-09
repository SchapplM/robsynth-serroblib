% Calculate vector of cutting forces with Newton-Euler
% S6PRPRPR2
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
% Datum: 2019-05-04 22:25
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRPR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 22:20:17
% EndTime: 2019-05-04 22:20:21
% DurationCPUTime: 2.76s
% Computational Cost: add. (40508->149), mult. (78880->206), div. (0->0), fcn. (55644->14), ass. (0->89)
t93 = cos(qJ(4));
t112 = t93 * qJD(2);
t84 = sin(pkin(6));
t94 = cos(qJ(2));
t116 = t84 * t94;
t83 = sin(pkin(10));
t87 = cos(pkin(10));
t70 = t83 * g(1) - t87 * g(2);
t88 = cos(pkin(6));
t118 = t70 * t88;
t71 = -t87 * g(1) - t83 * g(2);
t80 = -g(3) + qJDD(1);
t91 = sin(qJ(2));
t101 = t80 * t116 + t94 * t118 - t91 * t71;
t43 = qJDD(2) * pkin(2) + t101;
t117 = t84 * t91;
t108 = t80 * t117 + t91 * t118 + t94 * t71;
t96 = qJD(2) ^ 2;
t44 = -t96 * pkin(2) + t108;
t82 = sin(pkin(11));
t86 = cos(pkin(11));
t114 = t82 * t43 + t86 * t44;
t32 = -t96 * pkin(3) + qJDD(2) * pkin(8) + t114;
t103 = -t84 * t70 + t88 * t80;
t50 = qJDD(3) + t103;
t90 = sin(qJ(4));
t115 = t93 * t32 + t90 * t50;
t66 = (-pkin(4) * t93 - qJ(5) * t90) * qJD(2);
t95 = qJD(4) ^ 2;
t25 = -t95 * pkin(4) + qJDD(4) * qJ(5) + t66 * t112 + t115;
t111 = qJD(2) * qJD(4);
t106 = t93 * t111;
t104 = t86 * t43 - t82 * t44;
t31 = -qJDD(2) * pkin(3) - t96 * pkin(8) - t104;
t68 = t90 * qJDD(2) + t106;
t78 = t90 * t111;
t69 = t93 * qJDD(2) - t78;
t28 = (-t68 - t106) * qJ(5) + (-t69 + t78) * pkin(4) + t31;
t113 = qJD(2) * t90;
t81 = sin(pkin(12));
t85 = cos(pkin(12));
t63 = t81 * qJD(4) + t85 * t113;
t102 = -0.2e1 * qJD(5) * t63 - t81 * t25 + t85 * t28;
t54 = t81 * qJDD(4) + t85 * t68;
t62 = t85 * qJD(4) - t81 * t113;
t19 = (-t62 * t112 - t54) * pkin(9) + (t62 * t63 - t69) * pkin(5) + t102;
t110 = 0.2e1 * qJD(5) * t62 + t85 * t25 + t81 * t28;
t53 = t85 * qJDD(4) - t81 * t68;
t55 = -pkin(5) * t112 - t63 * pkin(9);
t61 = t62 ^ 2;
t20 = -t61 * pkin(5) + t53 * pkin(9) + t55 * t112 + t110;
t89 = sin(qJ(6));
t92 = cos(qJ(6));
t45 = t92 * t62 - t89 * t63;
t35 = t45 * qJD(6) + t89 * t53 + t92 * t54;
t46 = t89 * t62 + t92 * t63;
t37 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t76 = qJD(6) - t112;
t41 = -t76 * mrSges(7,2) + t45 * mrSges(7,3);
t64 = qJDD(6) - t69;
t15 = m(7) * (t92 * t19 - t89 * t20) - t35 * mrSges(7,3) + t64 * mrSges(7,1) - t46 * t37 + t76 * t41;
t34 = -t46 * qJD(6) + t92 * t53 - t89 * t54;
t42 = t76 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t89 * t19 + t92 * t20) + t34 * mrSges(7,3) - t64 * mrSges(7,2) + t45 * t37 - t76 * t42;
t47 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t51 = mrSges(6,2) * t112 + t62 * mrSges(6,3);
t13 = m(6) * t102 - t69 * mrSges(6,1) - t54 * mrSges(6,3) - t51 * t112 + t92 * t15 + t89 * t16 - t63 * t47;
t52 = -mrSges(6,1) * t112 - t63 * mrSges(6,3);
t14 = m(6) * t110 + t69 * mrSges(6,2) + t53 * mrSges(6,3) + t52 * t112 - t89 * t15 + t92 * t16 + t62 * t47;
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t113;
t73 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t112;
t119 = m(5) * t31 - t69 * mrSges(5,1) + t68 * mrSges(5,2) + t85 * t13 + t81 * t14 + (t90 * t72 - t93 * t73) * qJD(2);
t67 = (-mrSges(5,1) * t93 + mrSges(5,2) * t90) * qJD(2);
t12 = m(5) * t115 - qJDD(4) * mrSges(5,2) + t69 * mrSges(5,3) - qJD(4) * t72 + t67 * t112 - t81 * t13 + t85 * t14;
t105 = -t90 * t32 + t93 * t50;
t24 = -qJDD(4) * pkin(4) - t95 * qJ(5) + t66 * t113 + qJDD(5) - t105;
t99 = t34 * mrSges(7,1) + t45 * t41 - m(7) * (-t53 * pkin(5) - t61 * pkin(9) + t63 * t55 + t24) - t35 * mrSges(7,2) - t46 * t42;
t97 = m(6) * t24 - t53 * mrSges(6,1) + t54 * mrSges(6,2) - t62 * t51 + t63 * t52 - t99;
t18 = m(5) * t105 + qJDD(4) * mrSges(5,1) - t68 * mrSges(5,3) + qJD(4) * t73 - t67 * t113 - t97;
t109 = m(4) * t50 + t90 * t12 + t93 * t18;
t10 = m(4) * t104 + qJDD(2) * mrSges(4,1) - t96 * mrSges(4,2) - t119;
t7 = m(4) * t114 - t96 * mrSges(4,1) - qJDD(2) * mrSges(4,2) + t93 * t12 - t90 * t18;
t5 = m(3) * t101 + qJDD(2) * mrSges(3,1) - t96 * mrSges(3,2) + t86 * t10 + t82 * t7;
t6 = m(3) * t108 - t96 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t82 * t10 + t86 * t7;
t9 = m(3) * t103 + t109;
t107 = m(2) * t80 + t5 * t116 + t6 * t117 + t88 * t9;
t2 = m(2) * t71 - t91 * t5 + t94 * t6;
t1 = m(2) * t70 - t84 * t9 + (t5 * t94 + t6 * t91) * t88;
t3 = [-m(1) * g(1) - t83 * t1 + t87 * t2, t2, t6, t7, t12, t14, t16; -m(1) * g(2) + t87 * t1 + t83 * t2, t1, t5, t10, t18, t13, t15; -m(1) * g(3) + t107, t107, t9, t109, t119, t97, -t99;];
f_new  = t3;
