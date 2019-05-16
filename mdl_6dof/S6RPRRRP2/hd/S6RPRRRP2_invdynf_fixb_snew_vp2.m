% Calculate vector of cutting forces with Newton-Euler
% S6RPRRRP2
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-05-06 01:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRRRP2_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP2_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-06 01:12:33
% EndTime: 2019-05-06 01:12:38
% DurationCPUTime: 2.20s
% Computational Cost: add. (27972->174), mult. (54251->219), div. (0->0), fcn. (35528->10), ass. (0->87)
t89 = sin(qJ(3));
t117 = qJD(1) * t89;
t88 = sin(qJ(4));
t92 = cos(qJ(4));
t69 = t88 * qJD(3) + t117 * t92;
t115 = qJD(1) * qJD(3);
t93 = cos(qJ(3));
t107 = t93 * t115;
t74 = t89 * qJDD(1) + t107;
t51 = -t69 * qJD(4) + t92 * qJDD(3) - t88 * t74;
t68 = t92 * qJD(3) - t117 * t88;
t52 = t68 * qJD(4) + t88 * qJDD(3) + t92 * t74;
t87 = sin(qJ(5));
t91 = cos(qJ(5));
t55 = t87 * t68 + t91 * t69;
t27 = -t55 * qJD(5) + t91 * t51 - t87 * t52;
t54 = t91 * t68 - t87 * t69;
t28 = t54 * qJD(5) + t87 * t51 + t91 * t52;
t116 = t93 * qJD(1);
t79 = qJD(4) - t116;
t78 = qJD(5) + t79;
t43 = t78 * pkin(5) - t55 * qJ(6);
t44 = t78 * mrSges(7,1) - t55 * mrSges(7,3);
t53 = t54 ^ 2;
t90 = sin(qJ(1));
t94 = cos(qJ(1));
t110 = t90 * g(1) - t94 * g(2);
t70 = qJDD(1) * pkin(1) + t110;
t102 = -t94 * g(1) - t90 * g(2);
t96 = qJD(1) ^ 2;
t72 = -t96 * pkin(1) + t102;
t85 = sin(pkin(10));
t86 = cos(pkin(10));
t118 = t85 * t70 + t86 * t72;
t49 = -t96 * pkin(2) + qJDD(1) * pkin(7) + t118;
t84 = -g(3) + qJDD(2);
t104 = -t89 * t49 + t93 * t84;
t73 = (-pkin(3) * t93 - pkin(8) * t89) * qJD(1);
t95 = qJD(3) ^ 2;
t36 = -qJDD(3) * pkin(3) - t95 * pkin(8) + t73 * t117 - t104;
t59 = t79 * pkin(4) - t69 * pkin(9);
t66 = t68 ^ 2;
t98 = -t51 * pkin(4) - t66 * pkin(9) + t69 * t59 + t36;
t111 = m(7) * (-t27 * pkin(5) - t53 * qJ(6) + t55 * t43 + qJDD(6) + t98) + t28 * mrSges(7,2) + t55 * t44;
t41 = -t78 * mrSges(7,2) + t54 * mrSges(7,3);
t42 = -t78 * mrSges(6,2) + t54 * mrSges(6,3);
t45 = t78 * mrSges(6,1) - t55 * mrSges(6,3);
t127 = m(6) * t98 + t28 * mrSges(6,2) + t55 * t45 + t111 - (t42 + t41) * t54 - (mrSges(6,1) + mrSges(7,1)) * t27;
t57 = -t79 * mrSges(5,2) + t68 * mrSges(5,3);
t58 = t79 * mrSges(5,1) - t69 * mrSges(5,3);
t126 = m(5) * t36 - t51 * mrSges(5,1) + t52 * mrSges(5,2) - t68 * t57 + t69 * t58 + t127;
t103 = t86 * t70 - t85 * t72;
t48 = -qJDD(1) * pkin(2) - t96 * pkin(7) - t103;
t81 = t89 * t115;
t75 = t93 * qJDD(1) - t81;
t32 = (-t74 - t107) * pkin(8) + (-t75 + t81) * pkin(3) + t48;
t119 = t93 * t49 + t89 * t84;
t37 = -t95 * pkin(3) + qJDD(3) * pkin(8) + t116 * t73 + t119;
t105 = t92 * t32 - t88 * t37;
t67 = qJDD(4) - t75;
t19 = (t68 * t79 - t52) * pkin(9) + (t68 * t69 + t67) * pkin(4) + t105;
t121 = t88 * t32 + t92 * t37;
t21 = -t66 * pkin(4) + t51 * pkin(9) - t79 * t59 + t121;
t122 = t87 * t19 + t91 * t21;
t38 = -t54 * mrSges(7,1) + t55 * mrSges(7,2);
t112 = m(7) * (-t53 * pkin(5) + t27 * qJ(6) + 0.2e1 * qJD(6) * t54 - t78 * t43 + t122) + t27 * mrSges(7,3) + t54 * t38;
t39 = -t54 * mrSges(6,1) + t55 * mrSges(6,2);
t65 = qJDD(5) + t67;
t10 = m(6) * t122 + t27 * mrSges(6,3) + t54 * t39 + (-t45 - t44) * t78 + (-mrSges(6,2) - mrSges(7,2)) * t65 + t112;
t56 = -t68 * mrSges(5,1) + t69 * mrSges(5,2);
t106 = t91 * t19 - t87 * t21;
t113 = m(7) * (-0.2e1 * qJD(6) * t55 + (t54 * t78 - t28) * qJ(6) + (t54 * t55 + t65) * pkin(5) + t106) + t78 * t41 + t65 * mrSges(7,1);
t9 = m(6) * t106 + t65 * mrSges(6,1) + t78 * t42 + (-t39 - t38) * t55 + (-mrSges(6,3) - mrSges(7,3)) * t28 + t113;
t7 = m(5) * t105 + t67 * mrSges(5,1) - t52 * mrSges(5,3) + t87 * t10 - t69 * t56 + t79 * t57 + t91 * t9;
t76 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t117;
t77 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t116;
t8 = m(5) * t121 - t67 * mrSges(5,2) + t51 * mrSges(5,3) + t91 * t10 + t68 * t56 - t79 * t58 - t87 * t9;
t124 = m(4) * t48 - t75 * mrSges(4,1) + t74 * mrSges(4,2) + t92 * t7 + t88 * t8 + (t89 * t76 - t93 * t77) * qJD(1);
t71 = (-mrSges(4,1) * t93 + mrSges(4,2) * t89) * qJD(1);
t12 = m(4) * t104 + qJDD(3) * mrSges(4,1) - t74 * mrSges(4,3) + qJD(3) * t77 - t71 * t117 - t126;
t6 = m(4) * t119 - qJDD(3) * mrSges(4,2) + t75 * mrSges(4,3) - qJD(3) * t76 + t116 * t71 - t88 * t7 + t92 * t8;
t114 = m(3) * t84 + t93 * t12 + t89 * t6;
t4 = m(3) * t103 + qJDD(1) * mrSges(3,1) - t96 * mrSges(3,2) - t124;
t3 = m(3) * t118 - t96 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t89 * t12 + t93 * t6;
t2 = m(2) * t102 - t96 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t86 * t3 - t85 * t4;
t1 = m(2) * t110 + qJDD(1) * mrSges(2,1) - t96 * mrSges(2,2) + t85 * t3 + t86 * t4;
t5 = [-m(1) * g(1) - t90 * t1 + t94 * t2, t2, t3, t6, t8, t10, -t65 * mrSges(7,2) - t78 * t44 + t112; -m(1) * g(2) + t94 * t1 + t90 * t2, t1, t4, t12, t7, t9, -t28 * mrSges(7,3) - t55 * t38 + t113; (-m(1) - m(2)) * g(3) + t114, -m(2) * g(3) + t114, t114, t124, t126, t127, -t27 * mrSges(7,1) - t54 * t41 + t111;];
f_new  = t5;
