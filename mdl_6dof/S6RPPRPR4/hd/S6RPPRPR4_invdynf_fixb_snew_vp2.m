% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR4
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-05-05 14:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR4_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:14:55
% EndTime: 2019-05-05 14:14:58
% DurationCPUTime: 1.29s
% Computational Cost: add. (16567->155), mult. (32798->201), div. (0->0), fcn. (17939->10), ass. (0->81)
t68 = sin(pkin(10));
t70 = cos(pkin(10));
t73 = sin(qJ(4));
t76 = cos(qJ(4));
t51 = (t68 * t73 - t70 * t76) * qJD(1);
t99 = qJD(1) * qJD(4);
t96 = t76 * t99;
t54 = -t73 * qJDD(1) - t96;
t55 = -t76 * qJDD(1) + t73 * t99;
t101 = qJD(1) * t73;
t57 = qJD(4) * mrSges(5,1) + mrSges(5,3) * t101;
t100 = qJD(1) * t76;
t58 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t100;
t79 = qJD(1) ^ 2;
t52 = (t68 * t76 + t70 * t73) * qJD(1);
t33 = -t51 * pkin(5) + t52 * pkin(8);
t78 = qJD(4) ^ 2;
t107 = 2 * qJD(5);
t105 = -pkin(1) - pkin(2);
t74 = sin(qJ(1));
t77 = cos(qJ(1));
t92 = -t77 * g(1) - t74 * g(2);
t87 = qJDD(1) * qJ(2) + 0.2e1 * qJD(2) * qJD(1) + t92;
t43 = t105 * t79 + t87;
t97 = t74 * g(1) - t77 * g(2);
t84 = -t79 * qJ(2) + qJDD(2) - t97;
t47 = t105 * qJDD(1) + t84;
t69 = sin(pkin(9));
t71 = cos(pkin(9));
t102 = t71 * t43 + t69 * t47;
t27 = -t79 * pkin(3) - qJDD(1) * pkin(7) + t102;
t66 = g(3) + qJDD(3);
t95 = -t73 * t27 + t76 * t66;
t19 = (-t54 - t96) * qJ(5) + (t73 * t76 * t79 + qJDD(4)) * pkin(4) + t95;
t103 = t76 * t27 + t73 * t66;
t56 = qJD(4) * pkin(4) + qJ(5) * t101;
t65 = t76 ^ 2;
t20 = -t65 * t79 * pkin(4) + t55 * qJ(5) - qJD(4) * t56 + t103;
t98 = t51 * t107 + t68 * t19 + t70 * t20;
t15 = -t78 * pkin(5) + qJDD(4) * pkin(8) + t51 * t33 + t98;
t36 = -t68 * t54 + t70 * t55;
t37 = t70 * t54 + t68 * t55;
t94 = -t69 * t43 + t71 * t47;
t93 = qJDD(1) * pkin(3) - t94;
t80 = -t56 * t101 - t55 * pkin(4) + qJDD(5) + (-qJ(5) * t65 - pkin(7)) * t79 + t93;
t16 = (-qJD(4) * t51 - t37) * pkin(8) + (-qJD(4) * t52 - t36) * pkin(5) + t80;
t72 = sin(qJ(6));
t75 = cos(qJ(6));
t40 = t75 * qJD(4) + t72 * t52;
t24 = t40 * qJD(6) + t72 * qJDD(4) + t75 * t37;
t41 = t72 * qJD(4) - t75 * t52;
t28 = -t40 * mrSges(7,1) + t41 * mrSges(7,2);
t49 = qJD(6) - t51;
t29 = -t49 * mrSges(7,2) + t40 * mrSges(7,3);
t35 = qJDD(6) - t36;
t12 = m(7) * (-t72 * t15 + t75 * t16) - t24 * mrSges(7,3) + t35 * mrSges(7,1) - t41 * t28 + t49 * t29;
t23 = -t41 * qJD(6) + t75 * qJDD(4) - t72 * t37;
t30 = t49 * mrSges(7,1) - t41 * mrSges(7,3);
t13 = m(7) * (t75 * t15 + t72 * t16) + t23 * mrSges(7,3) - t35 * mrSges(7,2) + t40 * t28 - t49 * t30;
t44 = -qJD(4) * mrSges(6,2) + t51 * mrSges(6,3);
t45 = qJD(4) * mrSges(6,1) + t52 * mrSges(6,3);
t83 = -m(6) * t80 + t36 * mrSges(6,1) - t37 * mrSges(6,2) - t75 * t12 - t72 * t13 + t51 * t44 + t52 * t45;
t108 = -(t73 * t57 - t76 * t58) * qJD(1) + m(5) * (-t79 * pkin(7) + t93) - t55 * mrSges(5,1) + t54 * mrSges(5,2) - t83;
t106 = -m(2) - m(3);
t104 = mrSges(2,1) + mrSges(3,1);
t91 = -t70 * t19 + t68 * t20;
t53 = (mrSges(5,1) * t76 - mrSges(5,2) * t73) * qJD(1);
t32 = -t51 * mrSges(6,1) - t52 * mrSges(6,2);
t8 = m(6) * t98 - qJDD(4) * mrSges(6,2) + t36 * mrSges(6,3) - qJD(4) * t45 - t72 * t12 + t75 * t13 + t51 * t32;
t82 = m(7) * (-qJDD(4) * pkin(5) - t78 * pkin(8) + (-(2 * qJD(5)) - t33) * t52 + t91) - t23 * mrSges(7,1) + t24 * mrSges(7,2) - t40 * t29 + t41 * t30;
t9 = m(6) * (t52 * t107 - t91) - t37 * mrSges(6,3) + qJDD(4) * mrSges(6,1) + t52 * t32 + qJD(4) * t44 - t82;
t5 = m(5) * t95 + qJDD(4) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(4) * t58 + t53 * t101 + t68 * t8 + t70 * t9;
t6 = m(5) * t103 - qJDD(4) * mrSges(5,2) + t55 * mrSges(5,3) - qJD(4) * t57 - t53 * t100 - t68 * t9 + t70 * t8;
t4 = m(4) * t102 - t79 * mrSges(4,1) + qJDD(1) * mrSges(4,2) - t73 * t5 + t76 * t6;
t7 = m(4) * t94 - qJDD(1) * mrSges(4,1) - t79 * mrSges(4,2) - t108;
t88 = t71 * t4 - t69 * t7 + m(3) * (-t79 * pkin(1) + t87) + qJDD(1) * mrSges(3,3);
t86 = -m(3) * (-qJDD(1) * pkin(1) + t84) - t69 * t4 - t71 * t7;
t85 = m(4) * t66 + t76 * t5 + t73 * t6;
t2 = m(2) * t97 + (-mrSges(2,2) + mrSges(3,3)) * t79 + t104 * qJDD(1) + t86;
t1 = m(2) * t92 - qJDD(1) * mrSges(2,2) - t104 * t79 + t88;
t3 = [-m(1) * g(1) + t77 * t1 - t74 * t2, t1, -t79 * mrSges(3,1) + t88, t4, t6, t8, t13; -m(1) * g(2) + t74 * t1 + t77 * t2, t2, -m(3) * g(3) - t85, t7, t5, t9, t12; (-m(1) + t106) * g(3) - t85, t106 * g(3) - t85, -qJDD(1) * mrSges(3,1) - t79 * mrSges(3,3) - t86, t85, t108, -t83, t82;];
f_new  = t3;
