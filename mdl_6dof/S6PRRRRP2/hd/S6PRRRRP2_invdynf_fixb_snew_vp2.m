% Calculate vector of cutting forces with Newton-Euler
% S6PRRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-05-05 09:40
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 09:34:37
% EndTime: 2019-05-05 09:34:44
% DurationCPUTime: 2.54s
% Computational Cost: add. (34836->175), mult. (68087->229), div. (0->0), fcn. (48450->12), ass. (0->92)
t94 = sin(qJ(4));
t95 = sin(qJ(3));
t97 = cos(qJ(4));
t98 = cos(qJ(3));
t67 = (t94 * t95 - t97 * t98) * qJD(2);
t89 = sin(pkin(11));
t91 = cos(pkin(11));
t75 = t89 * g(1) - t91 * g(2);
t92 = cos(pkin(6));
t128 = t75 * t92;
t76 = -t91 * g(1) - t89 * g(2);
t88 = -g(3) + qJDD(1);
t90 = sin(pkin(6));
t96 = sin(qJ(2));
t99 = cos(qJ(2));
t134 = (t88 * t90 + t128) * t99 - t96 * t76;
t100 = qJD(2) ^ 2;
t103 = -qJDD(2) * pkin(2) - t134;
t119 = qJD(2) * t95;
t117 = qJD(2) * qJD(3);
t74 = t98 * qJDD(2) - t95 * t117;
t80 = qJD(3) * pkin(3) - pkin(9) * t119;
t87 = t98 ^ 2;
t101 = -t74 * pkin(3) + t80 * t119 + (-pkin(9) * t87 - pkin(8)) * t100 + t103;
t130 = cos(qJ(5));
t127 = t90 * t96;
t114 = t88 * t127 + t96 * t128 + t99 * t76;
t53 = -t100 * pkin(2) + qJDD(2) * pkin(8) + t114;
t63 = -t90 * t75 + t92 * t88;
t110 = -t95 * t53 + t98 * t63;
t112 = t98 * t117;
t73 = t95 * qJDD(2) + t112;
t29 = (-t73 + t112) * pkin(9) + (t100 * t95 * t98 + qJDD(3)) * pkin(3) + t110;
t121 = t98 * t53 + t95 * t63;
t30 = -t87 * t100 * pkin(3) + t74 * pkin(9) - qJD(3) * t80 + t121;
t123 = t94 * t29 + t97 * t30;
t68 = (t94 * t98 + t95 * t97) * qJD(2);
t56 = t67 * pkin(4) - t68 * pkin(10);
t86 = qJD(3) + qJD(4);
t84 = t86 ^ 2;
t85 = qJDD(3) + qJDD(4);
t23 = -t84 * pkin(4) + t85 * pkin(10) - t67 * t56 + t123;
t46 = -t68 * qJD(4) - t94 * t73 + t97 * t74;
t47 = -t67 * qJD(4) + t97 * t73 + t94 * t74;
t25 = (t67 * t86 - t47) * pkin(10) + (t68 * t86 - t46) * pkin(4) + t101;
t93 = sin(qJ(5));
t124 = t130 * t23 + t93 * t25;
t57 = -t130 * t86 + t93 * t68;
t58 = t130 * t68 + t93 * t86;
t38 = t57 * pkin(5) - t58 * qJ(6);
t44 = qJDD(5) - t46;
t65 = qJD(5) + t67;
t51 = -t65 * mrSges(7,1) + t58 * mrSges(7,2);
t64 = t65 ^ 2;
t116 = m(7) * (-t64 * pkin(5) + t44 * qJ(6) + 0.2e1 * qJD(6) * t65 - t57 * t38 + t124) + t65 * t51 + t44 * mrSges(7,3);
t39 = t57 * mrSges(7,1) - t58 * mrSges(7,3);
t122 = -t57 * mrSges(6,1) - t58 * mrSges(6,2) - t39;
t125 = -mrSges(6,3) - mrSges(7,2);
t32 = t58 * qJD(5) - t130 * t85 + t93 * t47;
t50 = t65 * mrSges(6,1) - t58 * mrSges(6,3);
t14 = m(6) * t124 - t44 * mrSges(6,2) + t122 * t57 + t125 * t32 - t65 * t50 + t116;
t106 = t130 * t25 - t93 * t23;
t131 = m(7) * (-t44 * pkin(5) - t64 * qJ(6) + t58 * t38 + qJDD(6) - t106);
t33 = -t57 * qJD(5) + t130 * t47 + t93 * t85;
t48 = -t57 * mrSges(7,2) + t65 * mrSges(7,3);
t49 = -t65 * mrSges(6,2) - t57 * mrSges(6,3);
t16 = m(6) * t106 - t131 + (t49 + t48) * t65 + t122 * t58 + (mrSges(6,1) + mrSges(7,1)) * t44 + t125 * t33;
t61 = -t86 * mrSges(5,2) - t67 * mrSges(5,3);
t62 = t86 * mrSges(5,1) - t68 * mrSges(5,3);
t105 = -m(5) * t101 + t46 * mrSges(5,1) - t47 * mrSges(5,2) - t130 * t16 - t93 * t14 - t67 * t61 - t68 * t62;
t77 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t119;
t118 = qJD(2) * t98;
t78 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t118;
t133 = (t95 * t77 - t98 * t78) * qJD(2) + m(4) * (-t100 * pkin(8) + t103) - t74 * mrSges(4,1) + t73 * mrSges(4,2) - t105;
t111 = t97 * t29 - t94 * t30;
t22 = -t85 * pkin(4) - t84 * pkin(10) + t68 * t56 - t111;
t115 = m(7) * (-0.2e1 * qJD(6) * t58 + (t57 * t65 - t33) * qJ(6) + (t58 * t65 + t32) * pkin(5) + t22) + t32 * mrSges(7,1) + t57 * t48;
t132 = m(6) * t22 + t32 * mrSges(6,1) + (t50 - t51) * t58 + (mrSges(6,2) - mrSges(7,3)) * t33 + t57 * t49 + t115;
t10 = m(3) * t134 + qJDD(2) * mrSges(3,1) - t100 * mrSges(3,2) - t133;
t129 = t10 * t99;
t55 = t67 * mrSges(5,1) + t68 * mrSges(5,2);
t11 = m(5) * t123 - t85 * mrSges(5,2) + t46 * mrSges(5,3) + t130 * t14 - t93 * t16 - t67 * t55 - t86 * t62;
t12 = m(5) * t111 + t85 * mrSges(5,1) - t47 * mrSges(5,3) - t68 * t55 + t86 * t61 - t132;
t72 = (-mrSges(4,1) * t98 + mrSges(4,2) * t95) * qJD(2);
t7 = m(4) * t110 + qJDD(3) * mrSges(4,1) - t73 * mrSges(4,3) + qJD(3) * t78 + t94 * t11 - t72 * t119 + t97 * t12;
t8 = m(4) * t121 - qJDD(3) * mrSges(4,2) + t74 * mrSges(4,3) - qJD(3) * t77 + t97 * t11 + t72 * t118 - t94 * t12;
t4 = m(3) * t114 - t100 * mrSges(3,1) - qJDD(2) * mrSges(3,2) - t95 * t7 + t98 * t8;
t6 = m(3) * t63 + t98 * t7 + t95 * t8;
t113 = m(2) * t88 + t4 * t127 + t90 * t129 + t92 * t6;
t2 = m(2) * t76 - t96 * t10 + t99 * t4;
t1 = m(2) * t75 - t90 * t6 + (t4 * t96 + t129) * t92;
t3 = [-m(1) * g(1) - t89 * t1 + t91 * t2, t2, t4, t8, t11, t14, -t32 * mrSges(7,2) - t57 * t39 + t116; -m(1) * g(2) + t91 * t1 + t89 * t2, t1, t10, t7, t12, t16, -t33 * mrSges(7,3) - t58 * t51 + t115; -m(1) * g(3) + t113, t113, t6, t133, -t105, t132, -t44 * mrSges(7,1) + t33 * mrSges(7,2) + t58 * t39 - t65 * t48 + t131;];
f_new  = t3;
