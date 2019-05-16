% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR2
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
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2,theta3]';
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
% Datum: 2019-05-05 15:22
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:20:22
% EndTime: 2019-05-05 15:20:27
% DurationCPUTime: 3.19s
% Computational Cost: add. (41646->169), mult. (93198->219), div. (0->0), fcn. (66740->12), ass. (0->95)
t101 = qJD(1) ^ 2;
t90 = cos(pkin(11));
t86 = t90 ^ 2;
t88 = sin(pkin(11));
t121 = t88 ^ 2 + t86;
t128 = mrSges(4,3) * t121;
t127 = pkin(3) * t101;
t94 = sin(qJ(4));
t126 = t88 * t94;
t100 = qJD(4) ^ 2;
t118 = qJD(1) * qJD(3);
t87 = -g(3) + qJDD(2);
t122 = -0.2e1 * t88 * t118 + t90 * t87;
t95 = sin(qJ(1));
t99 = cos(qJ(1));
t115 = t95 * g(1) - t99 * g(2);
t75 = qJDD(1) * pkin(1) + t115;
t111 = -t99 * g(1) - t95 * g(2);
t76 = -t101 * pkin(1) + t111;
t89 = sin(pkin(10));
t91 = cos(pkin(10));
t123 = t89 * t75 + t91 * t76;
t52 = -t101 * pkin(2) + qJDD(1) * qJ(3) + t123;
t37 = (-pkin(7) * qJDD(1) + t127 * t90 - t52) * t88 + t122;
t116 = t88 * t87 + (0.2e1 * t118 + t52) * t90;
t119 = qJDD(1) * t90;
t40 = pkin(7) * t119 - t127 * t86 + t116;
t98 = cos(qJ(4));
t124 = t94 * t37 + t98 * t40;
t69 = (t90 * t98 - t126) * qJD(1);
t108 = t88 * t98 + t90 * t94;
t70 = t108 * qJD(1);
t58 = -t69 * pkin(4) - t70 * pkin(8);
t26 = -t100 * pkin(4) + qJDD(4) * pkin(8) + t69 * t58 + t124;
t112 = t91 * t75 - t89 * t76;
t110 = qJDD(3) - t112;
t103 = (-pkin(3) * t90 - pkin(2)) * qJDD(1) + (-pkin(7) * t121 - qJ(3)) * t101 + t110;
t120 = t69 * qJD(4);
t67 = t70 * qJD(4);
t59 = -qJDD(1) * t126 + t119 * t98 - t67;
t60 = qJDD(1) * t108 + t120;
t29 = (-t60 - t120) * pkin(8) + (-t59 + t67) * pkin(4) + t103;
t93 = sin(qJ(5));
t97 = cos(qJ(5));
t125 = t97 * t26 + t93 * t29;
t109 = -t90 * mrSges(4,1) + t88 * mrSges(4,2);
t107 = qJDD(1) * mrSges(4,3) + t101 * t109;
t63 = t93 * qJD(4) + t97 * t70;
t38 = -t63 * qJD(5) + t97 * qJDD(4) - t93 * t60;
t62 = t97 * qJD(4) - t93 * t70;
t39 = t62 * qJD(5) + t93 * qJDD(4) + t97 * t60;
t92 = sin(qJ(6));
t96 = cos(qJ(6));
t44 = t92 * t62 + t96 * t63;
t23 = -t44 * qJD(6) + t96 * t38 - t92 * t39;
t43 = t96 * t62 - t92 * t63;
t24 = t43 * qJD(6) + t92 * t38 + t96 * t39;
t113 = t98 * t37 - t94 * t40;
t25 = -qJDD(4) * pkin(4) - t100 * pkin(8) + t70 * t58 - t113;
t68 = qJD(5) - t69;
t66 = qJD(6) + t68;
t32 = -t66 * mrSges(7,2) + t43 * mrSges(7,3);
t33 = t66 * mrSges(7,1) - t44 * mrSges(7,3);
t49 = t68 * pkin(5) - t63 * pkin(9);
t61 = t62 ^ 2;
t106 = t23 * mrSges(7,1) + t43 * t32 - m(7) * (-t38 * pkin(5) - t61 * pkin(9) + t63 * t49 + t25) - t24 * mrSges(7,2) - t44 * t33;
t47 = -t68 * mrSges(6,2) + t62 * mrSges(6,3);
t48 = t68 * mrSges(6,1) - t63 * mrSges(6,3);
t102 = m(6) * t25 - t38 * mrSges(6,1) + t39 * mrSges(6,2) - t62 * t47 + t63 * t48 - t106;
t55 = -t69 * mrSges(5,1) + t70 * mrSges(5,2);
t64 = -qJD(4) * mrSges(5,2) + t69 * mrSges(5,3);
t14 = m(5) * t113 + qJDD(4) * mrSges(5,1) - t60 * mrSges(5,3) + qJD(4) * t64 - t70 * t55 - t102;
t114 = -t93 * t26 + t97 * t29;
t57 = qJDD(5) - t59;
t17 = (t62 * t68 - t39) * pkin(9) + (t62 * t63 + t57) * pkin(5) + t114;
t18 = -t61 * pkin(5) + t38 * pkin(9) - t68 * t49 + t125;
t31 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t54 = qJDD(6) + t57;
t15 = m(7) * (t96 * t17 - t92 * t18) - t24 * mrSges(7,3) + t54 * mrSges(7,1) - t44 * t31 + t66 * t32;
t16 = m(7) * (t92 * t17 + t96 * t18) + t23 * mrSges(7,3) - t54 * mrSges(7,2) + t43 * t31 - t66 * t33;
t45 = -t62 * mrSges(6,1) + t63 * mrSges(6,2);
t12 = m(6) * t114 + t57 * mrSges(6,1) - t39 * mrSges(6,3) + t96 * t15 + t92 * t16 - t63 * t45 + t68 * t47;
t13 = m(6) * t125 - t57 * mrSges(6,2) + t38 * mrSges(6,3) - t92 * t15 + t96 * t16 + t62 * t45 - t68 * t48;
t65 = qJD(4) * mrSges(5,1) - t70 * mrSges(5,3);
t9 = m(5) * t124 - qJDD(4) * mrSges(5,2) + t59 * mrSges(5,3) - qJD(4) * t65 - t93 * t12 + t97 * t13 + t69 * t55;
t6 = m(4) * t122 + t94 * t9 + t98 * t14 + (-m(4) * t52 - t107) * t88;
t7 = m(4) * t116 + t107 * t90 - t94 * t14 + t98 * t9;
t117 = m(3) * t87 + t90 * t6 + t88 * t7;
t105 = -m(5) * t103 + t59 * mrSges(5,1) - t60 * mrSges(5,2) - t97 * t12 - t93 * t13 + t69 * t64 - t70 * t65;
t104 = m(4) * (-qJDD(1) * pkin(2) - t101 * qJ(3) + t110) - t105;
t8 = m(3) * t112 + (-mrSges(3,2) + t128) * t101 + (mrSges(3,1) - t109) * qJDD(1) - t104;
t3 = m(3) * t123 - t101 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t88 * t6 + t90 * t7;
t2 = m(2) * t111 - t101 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t91 * t3 - t89 * t8;
t1 = m(2) * t115 + qJDD(1) * mrSges(2,1) - t101 * mrSges(2,2) + t89 * t3 + t91 * t8;
t4 = [-m(1) * g(1) - t95 * t1 + t99 * t2, t2, t3, t7, t9, t13, t16; -m(1) * g(2) + t99 * t1 + t95 * t2, t1, t8, t6, t14, t12, t15; (-m(1) - m(2)) * g(3) + t117, -m(2) * g(3) + t117, t117, qJDD(1) * t109 - t101 * t128 + t104, -t105, t102, -t106;];
f_new  = t4;
