% Calculate vector of cutting forces with Newton-Euler
% S6RPPPRR2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 13:37
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPPRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 13:35:50
% EndTime: 2019-05-05 13:35:52
% DurationCPUTime: 1.10s
% Computational Cost: add. (12407->139), mult. (25813->176), div. (0->0), fcn. (16398->10), ass. (0->81)
t111 = -pkin(2) - qJ(4);
t80 = qJD(1) ^ 2;
t75 = sin(qJ(1));
t78 = cos(qJ(1));
t98 = t75 * g(1) - t78 * g(2);
t53 = qJDD(1) * pkin(1) + t98;
t94 = -t78 * g(1) - t75 * g(2);
t54 = -t80 * pkin(1) + t94;
t70 = sin(pkin(9));
t72 = cos(pkin(9));
t96 = t72 * t53 - t70 * t54;
t86 = -t80 * qJ(3) + qJDD(3) - t96;
t116 = -(2 * qJD(1) * qJD(4)) + t111 * qJDD(1) + t86;
t69 = sin(pkin(10));
t65 = t69 ^ 2;
t71 = cos(pkin(10));
t107 = t71 ^ 2 + t65;
t97 = t107 * mrSges(5,3);
t74 = sin(qJ(5));
t77 = cos(qJ(5));
t92 = t69 * t77 + t71 * t74;
t47 = t92 * qJD(1);
t108 = t70 * t53 + t72 * t54;
t115 = qJDD(1) * qJ(3) + (2 * qJD(3) * qJD(1)) + t108;
t91 = -t69 * t74 + t71 * t77;
t48 = t91 * qJD(1);
t105 = t48 * qJD(5);
t40 = -t92 * qJDD(1) - t105;
t114 = pkin(4) * t80;
t113 = mrSges(3,1) - mrSges(4,2);
t112 = -mrSges(3,2) + mrSges(4,3);
t103 = qJDD(1) * t71;
t109 = t116 * t71;
t68 = -g(3) + qJDD(2);
t19 = -pkin(7) * t103 + (-t71 * t114 - t68) * t69 + t109;
t100 = t116 * t69 + t71 * t68;
t104 = qJDD(1) * t69;
t20 = -pkin(7) * t104 - t65 * t114 + t100;
t110 = t74 * t19 + t77 * t20;
t106 = t47 * qJD(5);
t39 = t47 * pkin(5) - t48 * pkin(8);
t79 = qJD(5) ^ 2;
t15 = -t79 * pkin(5) + qJDD(5) * pkin(8) - t47 * t39 + t110;
t41 = t91 * qJDD(1) - t106;
t89 = qJDD(4) + t115;
t82 = pkin(4) * t104 + (-t107 * pkin(7) + t111) * t80 + t89;
t16 = (-t41 + t106) * pkin(8) + (-t40 + t105) * pkin(5) + t82;
t73 = sin(qJ(6));
t76 = cos(qJ(6));
t42 = t76 * qJD(5) - t73 * t48;
t22 = t42 * qJD(6) + t73 * qJDD(5) + t76 * t41;
t43 = t73 * qJD(5) + t76 * t48;
t25 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t46 = qJD(6) + t47;
t30 = -t46 * mrSges(7,2) + t42 * mrSges(7,3);
t38 = qJDD(6) - t40;
t12 = m(7) * (-t73 * t15 + t76 * t16) - t22 * mrSges(7,3) + t38 * mrSges(7,1) - t43 * t25 + t46 * t30;
t21 = -t43 * qJD(6) + t76 * qJDD(5) - t73 * t41;
t31 = t46 * mrSges(7,1) - t43 * mrSges(7,3);
t13 = m(7) * (t76 * t15 + t73 * t16) + t21 * mrSges(7,3) - t38 * mrSges(7,2) + t42 * t25 - t46 * t31;
t36 = t47 * mrSges(6,1) + t48 * mrSges(6,2);
t45 = qJD(5) * mrSges(6,1) - t48 * mrSges(6,3);
t8 = m(6) * t110 - qJDD(5) * mrSges(6,2) + t40 * mrSges(6,3) - qJD(5) * t45 - t73 * t12 + t76 * t13 - t47 * t36;
t88 = -qJDD(1) * mrSges(5,3) - t80 * (mrSges(5,1) * t69 + mrSges(5,2) * t71);
t44 = -qJD(5) * mrSges(6,2) - t47 * mrSges(6,3);
t93 = t77 * t19 - t74 * t20;
t84 = m(7) * (-qJDD(5) * pkin(5) - t79 * pkin(8) + t48 * t39 - t93) - t21 * mrSges(7,1) + t22 * mrSges(7,2) - t42 * t30 + t43 * t31;
t9 = m(6) * t93 + qJDD(5) * mrSges(6,1) - t41 * mrSges(6,3) + qJD(5) * t44 - t48 * t36 - t84;
t5 = m(5) * (-t69 * t68 + t109) + t74 * t8 + t77 * t9 + t88 * t71;
t6 = m(5) * t100 + t88 * t69 - t74 * t9 + t77 * t8;
t95 = m(4) * t68 - t69 * t5 + t71 * t6;
t90 = m(3) * t68 + t95;
t87 = -m(4) * (-qJDD(1) * pkin(2) + t86) - t71 * t5 - t69 * t6;
t85 = -m(6) * t82 + t40 * mrSges(6,1) - t41 * mrSges(6,2) - t76 * t12 - t73 * t13 - t47 * t44 - t48 * t45;
t83 = m(5) * (t111 * t80 + t89) + mrSges(5,1) * t104 + mrSges(5,2) * t103 - t85;
t81 = m(4) * (t80 * pkin(2) - t115) - t83;
t7 = -t81 + (-t97 - t113) * t80 + t112 * qJDD(1) + m(3) * t108;
t3 = m(3) * t96 + t113 * qJDD(1) + t112 * t80 + t87;
t2 = m(2) * t94 - t80 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t70 * t3 + t72 * t7;
t1 = m(2) * t98 + qJDD(1) * mrSges(2,1) - t80 * mrSges(2,2) + t72 * t3 + t70 * t7;
t4 = [-m(1) * g(1) - t75 * t1 + t78 * t2, t2, t7, t95, t6, t8, t13; -m(1) * g(2) + t78 * t1 + t75 * t2, t1, t3, t81 - qJDD(1) * mrSges(4,3) + (-mrSges(4,2) + t97) * t80, t5, t9, t12; (-m(1) - m(2)) * g(3) + t90, -m(2) * g(3) + t90, t90, qJDD(1) * mrSges(4,2) - t80 * mrSges(4,3) - t87, -t80 * t97 + t83, -t85, t84;];
f_new  = t4;
