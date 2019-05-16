% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR3
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-05-05 14:10
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR3_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:08:58
% EndTime: 2019-05-05 14:09:00
% DurationCPUTime: 1.19s
% Computational Cost: add. (13632->152), mult. (27691->198), div. (0->0), fcn. (16712->10), ass. (0->81)
t74 = sin(pkin(10));
t76 = cos(pkin(10));
t79 = sin(qJ(4));
t82 = cos(qJ(4));
t48 = (t74 * t82 + t76 * t79) * qJD(1);
t80 = sin(qJ(1));
t83 = cos(qJ(1));
t100 = t80 * g(1) - t83 * g(2);
t57 = qJDD(1) * pkin(1) + t100;
t85 = qJD(1) ^ 2;
t96 = -t83 * g(1) - t80 * g(2);
t59 = -t85 * pkin(1) + t96;
t75 = sin(pkin(9));
t77 = cos(pkin(9));
t108 = t75 * t57 + t77 * t59;
t98 = qJDD(1) * qJ(3) + 0.2e1 * qJD(3) * qJD(1) + t108;
t113 = 2 * qJD(5);
t112 = -pkin(2) - pkin(7);
t111 = pkin(4) * t85;
t110 = mrSges(3,1) - mrSges(4,2);
t109 = -mrSges(3,2) + mrSges(4,3);
t99 = t77 * t57 - t75 * t59;
t91 = -t85 * qJ(3) + qJDD(3) - t99;
t32 = t112 * qJDD(1) + t91;
t73 = -g(3) + qJDD(2);
t107 = t79 * t32 + t82 * t73;
t106 = qJD(1) * t79;
t105 = qJD(1) * t82;
t104 = qJD(1) * qJD(4);
t27 = t82 * t32;
t61 = t82 * qJDD(1) - t79 * t104;
t19 = qJDD(4) * pkin(4) - t61 * qJ(5) + t27 + (-qJ(5) * t104 - t82 * t111 - t73) * t79;
t60 = -t79 * qJDD(1) - t82 * t104;
t63 = qJD(4) * pkin(4) - qJ(5) * t105;
t72 = t79 ^ 2;
t20 = t60 * qJ(5) - qJD(4) * t63 - t72 * t111 + t107;
t102 = -t48 * t113 + t74 * t19 + t76 * t20;
t58 = (mrSges(5,1) * t79 + mrSges(5,2) * t82) * qJD(1);
t62 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t106;
t49 = t76 * t105 - t74 * t106;
t37 = t48 * pkin(5) - t49 * pkin(8);
t84 = qJD(4) ^ 2;
t15 = -t84 * pkin(5) + qJDD(4) * pkin(8) - t48 * t37 + t102;
t40 = t76 * t60 - t74 * t61;
t41 = t74 * t60 + t76 * t61;
t88 = -t60 * pkin(4) + qJDD(5) + t63 * t105 + (-qJ(5) * t72 + t112) * t85 + t98;
t16 = (qJD(4) * t48 - t41) * pkin(8) + (qJD(4) * t49 - t40) * pkin(5) + t88;
t78 = sin(qJ(6));
t81 = cos(qJ(6));
t42 = t81 * qJD(4) - t78 * t49;
t24 = t42 * qJD(6) + t78 * qJDD(4) + t81 * t41;
t43 = t78 * qJD(4) + t81 * t49;
t25 = -t42 * mrSges(7,1) + t43 * mrSges(7,2);
t47 = qJD(6) + t48;
t29 = -t47 * mrSges(7,2) + t42 * mrSges(7,3);
t39 = qJDD(6) - t40;
t12 = m(7) * (-t78 * t15 + t81 * t16) - t24 * mrSges(7,3) + t39 * mrSges(7,1) - t43 * t25 + t47 * t29;
t23 = -t43 * qJD(6) + t81 * qJDD(4) - t78 * t41;
t30 = t47 * mrSges(7,1) - t43 * mrSges(7,3);
t13 = m(7) * (t81 * t15 + t78 * t16) + t23 * mrSges(7,3) - t39 * mrSges(7,2) + t42 * t25 - t47 * t30;
t36 = t48 * mrSges(6,1) + t49 * mrSges(6,2);
t45 = qJD(4) * mrSges(6,1) - t49 * mrSges(6,3);
t8 = m(6) * t102 - qJDD(4) * mrSges(6,2) + t40 * mrSges(6,3) - qJD(4) * t45 - t78 * t12 + t81 * t13 - t48 * t36;
t44 = -qJD(4) * mrSges(6,2) - t48 * mrSges(6,3);
t95 = -t76 * t19 + t74 * t20;
t89 = m(7) * (-qJDD(4) * pkin(5) - t84 * pkin(8) + (t113 + t37) * t49 + t95) - t23 * mrSges(7,1) + t24 * mrSges(7,2) - t42 * t29 + t43 * t30;
t9 = m(6) * (-0.2e1 * qJD(5) * t49 - t95) - t41 * mrSges(6,3) + qJDD(4) * mrSges(6,1) - t49 * t36 + qJD(4) * t44 - t89;
t5 = m(5) * (-t79 * t73 + t27) - t61 * mrSges(5,3) + qJDD(4) * mrSges(5,1) - t58 * t105 + qJD(4) * t62 + t74 * t8 + t76 * t9;
t64 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t105;
t6 = m(5) * t107 - qJDD(4) * mrSges(5,2) + t60 * mrSges(5,3) - qJD(4) * t64 - t58 * t106 - t74 * t9 + t76 * t8;
t97 = m(4) * t73 - t79 * t5 + t82 * t6;
t93 = m(3) * t73 + t97;
t92 = -m(4) * (-qJDD(1) * pkin(2) + t91) - t82 * t5 - t79 * t6;
t90 = m(6) * t88 - t40 * mrSges(6,1) + t41 * mrSges(6,2) + t81 * t12 + t78 * t13 + t48 * t44 + t49 * t45;
t87 = -t60 * mrSges(5,1) + m(5) * (t112 * t85 + t98) + t62 * t106 + t64 * t105 + t61 * mrSges(5,2) + t90;
t86 = -m(4) * (t85 * pkin(2) - t98) + t87;
t7 = m(3) * t108 + t109 * qJDD(1) - t110 * t85 + t86;
t3 = m(3) * t99 + t110 * qJDD(1) + t109 * t85 + t92;
t2 = m(2) * t96 - t85 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t75 * t3 + t77 * t7;
t1 = m(2) * t100 + qJDD(1) * mrSges(2,1) - t85 * mrSges(2,2) + t77 * t3 + t75 * t7;
t4 = [-m(1) * g(1) - t80 * t1 + t83 * t2, t2, t7, t97, t6, t8, t13; -m(1) * g(2) + t83 * t1 + t80 * t2, t1, t3, -t85 * mrSges(4,2) - qJDD(1) * mrSges(4,3) - t86, t5, t9, t12; (-m(1) - m(2)) * g(3) + t93, -m(2) * g(3) + t93, t93, qJDD(1) * mrSges(4,2) - t85 * mrSges(4,3) - t92, t87, t90, t89;];
f_new  = t4;
