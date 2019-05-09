% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR1
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
% Datum: 2019-05-05 15:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 15:09:50
% EndTime: 2019-05-05 15:09:57
% DurationCPUTime: 3.33s
% Computational Cost: add. (44502->168), mult. (100949->218), div. (0->0), fcn. (74353->12), ass. (0->94)
t97 = qJD(1) ^ 2;
t87 = cos(pkin(11));
t82 = t87 ^ 2;
t85 = sin(pkin(11));
t119 = t85 ^ 2 + t82;
t125 = t119 * mrSges(4,3);
t124 = pkin(3) * t97;
t118 = pkin(7) * qJDD(1);
t116 = qJD(1) * qJD(3);
t84 = -g(3) + qJDD(2);
t120 = -0.2e1 * t85 * t116 + t87 * t84;
t92 = sin(qJ(1));
t96 = cos(qJ(1));
t113 = t92 * g(1) - t96 * g(2);
t71 = qJDD(1) * pkin(1) + t113;
t110 = -t96 * g(1) - t92 * g(2);
t72 = -t97 * pkin(1) + t110;
t86 = sin(pkin(10));
t88 = cos(pkin(10));
t121 = t86 * t71 + t88 * t72;
t54 = -t97 * pkin(2) + qJDD(1) * qJ(3) + t121;
t41 = (t87 * t124 - t118 - t54) * t85 + t120;
t114 = t85 * t84 + (0.2e1 * t116 + t54) * t87;
t42 = t87 * t118 - t82 * t124 + t114;
t91 = sin(qJ(4));
t95 = cos(qJ(4));
t112 = t95 * t41 - t91 * t42;
t105 = -t85 * t91 + t87 * t95;
t65 = t105 * qJD(1);
t117 = t65 * qJD(4);
t106 = t85 * t95 + t87 * t91;
t60 = t106 * qJDD(1) + t117;
t66 = t106 * qJD(1);
t21 = (-t60 + t117) * pkin(8) + (t65 * t66 + qJDD(4)) * pkin(4) + t112;
t122 = t91 * t41 + t95 * t42;
t59 = -t66 * qJD(4) + t105 * qJDD(1);
t63 = qJD(4) * pkin(4) - t66 * pkin(8);
t64 = t65 ^ 2;
t23 = -t64 * pkin(4) + t59 * pkin(8) - qJD(4) * t63 + t122;
t90 = sin(qJ(5));
t94 = cos(qJ(5));
t123 = t90 * t21 + t94 * t23;
t108 = -t87 * mrSges(4,1) + t85 * mrSges(4,2);
t104 = qJDD(1) * mrSges(4,3) + t97 * t108;
t52 = t94 * t65 - t90 * t66;
t53 = t90 * t65 + t94 * t66;
t37 = -t52 * pkin(5) - t53 * pkin(9);
t83 = qJD(4) + qJD(5);
t79 = t83 ^ 2;
t80 = qJDD(4) + qJDD(5);
t18 = -t79 * pkin(5) + t80 * pkin(9) + t52 * t37 + t123;
t30 = -t53 * qJD(5) + t94 * t59 - t90 * t60;
t31 = t52 * qJD(5) + t90 * t59 + t94 * t60;
t111 = t88 * t71 - t86 * t72;
t109 = qJDD(3) - t111;
t100 = (-pkin(3) * t87 - pkin(2)) * qJDD(1) + (-t119 * pkin(7) - qJ(3)) * t97 + t109;
t98 = -t59 * pkin(4) - t64 * pkin(8) + t66 * t63 + t100;
t19 = (-t52 * t83 - t31) * pkin(9) + (t53 * t83 - t30) * pkin(5) + t98;
t89 = sin(qJ(6));
t93 = cos(qJ(6));
t45 = -t89 * t53 + t93 * t83;
t25 = t45 * qJD(6) + t93 * t31 + t89 * t80;
t29 = qJDD(6) - t30;
t46 = t93 * t53 + t89 * t83;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t50 = qJD(6) - t52;
t33 = -t50 * mrSges(7,2) + t45 * mrSges(7,3);
t15 = m(7) * (-t89 * t18 + t93 * t19) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t46 * t32 + t50 * t33;
t24 = -t46 * qJD(6) - t89 * t31 + t93 * t80;
t34 = t50 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t93 * t18 + t89 * t19) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t45 * t32 - t50 * t34;
t36 = -t52 * mrSges(6,1) + t53 * mrSges(6,2);
t48 = t83 * mrSges(6,1) - t53 * mrSges(6,3);
t11 = m(6) * t123 - t80 * mrSges(6,2) + t30 * mrSges(6,3) - t89 * t15 + t93 * t16 + t52 * t36 - t83 * t48;
t107 = t94 * t21 - t90 * t23;
t102 = m(7) * (-t80 * pkin(5) - t79 * pkin(9) + t53 * t37 - t107) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t45 * t33 + t46 * t34;
t47 = -t83 * mrSges(6,2) + t52 * mrSges(6,3);
t12 = m(6) * t107 + t80 * mrSges(6,1) - t31 * mrSges(6,3) - t53 * t36 + t83 * t47 - t102;
t57 = -t65 * mrSges(5,1) + t66 * mrSges(5,2);
t61 = -qJD(4) * mrSges(5,2) + t65 * mrSges(5,3);
t8 = m(5) * t112 + qJDD(4) * mrSges(5,1) - t60 * mrSges(5,3) + qJD(4) * t61 + t90 * t11 + t94 * t12 - t66 * t57;
t62 = qJD(4) * mrSges(5,1) - t66 * mrSges(5,3);
t9 = m(5) * t122 - qJDD(4) * mrSges(5,2) + t59 * mrSges(5,3) - qJD(4) * t62 + t94 * t11 - t90 * t12 + t65 * t57;
t6 = m(4) * t120 + t91 * t9 + t95 * t8 + (-m(4) * t54 - t104) * t85;
t7 = m(4) * t114 + t104 * t87 - t91 * t8 + t95 * t9;
t115 = m(3) * t84 + t87 * t6 + t85 * t7;
t103 = -m(6) * t98 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t93 * t15 - t89 * t16 + t52 * t47 - t53 * t48;
t101 = -m(5) * t100 + t59 * mrSges(5,1) - t60 * mrSges(5,2) + t65 * t61 - t66 * t62 + t103;
t99 = m(4) * (-qJDD(1) * pkin(2) - t97 * qJ(3) + t109) - t101;
t10 = (-mrSges(3,2) + t125) * t97 + (mrSges(3,1) - t108) * qJDD(1) + m(3) * t111 - t99;
t3 = m(3) * t121 - t97 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t85 * t6 + t87 * t7;
t2 = m(2) * t110 - t97 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t86 * t10 + t88 * t3;
t1 = m(2) * t113 + qJDD(1) * mrSges(2,1) - t97 * mrSges(2,2) + t88 * t10 + t86 * t3;
t4 = [-m(1) * g(1) - t92 * t1 + t96 * t2, t2, t3, t7, t9, t11, t16; -m(1) * g(2) + t96 * t1 + t92 * t2, t1, t10, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t115, t108 * qJDD(1) - t97 * t125 + t99, -t101, -t103, t102;];
f_new  = t4;
