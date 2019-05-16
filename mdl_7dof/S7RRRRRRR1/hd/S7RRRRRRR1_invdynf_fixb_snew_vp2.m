% Calculate vector of cutting forces with Newton-Euler
% S7RRRRRRR1
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [7x1]
%   Generalized joint coordinates (joint angles)
% qJD [7x1]
%   Generalized joint velocities
% qJDD [7x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d1,d3,d5,d7]';
% m_mdh [8x1]
%   mass of all robot links (including the base)
% mrSges [8x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [8x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x8]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-05-09 00:12
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S7RRRRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(7,1),zeros(7,1),zeros(7,1),zeros(3,1),zeros(4,1),zeros(8,1),zeros(8,3),zeros(8,6)}
assert(isreal(qJ) && all(size(qJ) == [7 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: qJ has to be [7x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [7 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: qJD has to be [7x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [7 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [7x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: pkin has to be [4x1] (double)');
assert(isreal(m) && all(size(m) == [8 1]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: m has to be [8x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [8,3]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [8x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [8 6]), ...
  'S7RRRRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [8x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-08 21:40:28
% EndTime: 2019-05-08 21:40:55
% DurationCPUTime: 8.09s
% Computational Cost: add. (152986->206), mult. (276728->274), div. (0->0), fcn. (230115->14), ass. (0->111)
t105 = sin(qJ(3));
t106 = sin(qJ(2));
t112 = cos(qJ(3));
t113 = cos(qJ(2));
t102 = sin(qJ(6));
t109 = cos(qJ(6));
t101 = sin(qJ(7));
t108 = cos(qJ(7));
t103 = sin(qJ(5));
t110 = cos(qJ(5));
t104 = sin(qJ(4));
t111 = cos(qJ(4));
t131 = qJD(1) * t106;
t88 = -t112 * qJD(2) - t105 * t131;
t129 = qJD(1) * qJD(2);
t127 = t113 * t129;
t92 = t106 * qJDD(1) + t127;
t74 = t88 * qJD(3) - t105 * qJDD(2) + t112 * t92;
t89 = -t105 * qJD(2) + t112 * t131;
t130 = t113 * qJD(1);
t98 = qJD(3) + t130;
t81 = -t104 * t89 - t111 * t98;
t93 = t113 * qJDD(1) - t106 * t129;
t87 = qJDD(3) + t93;
t54 = t81 * qJD(4) - t104 * t87 + t111 * t74;
t107 = sin(qJ(1));
t114 = cos(qJ(1));
t96 = t107 * g(1) - t114 * g(2);
t77 = (t92 + t127) * pkin(2) - t96;
t115 = qJD(1) ^ 2;
t97 = -t114 * g(1) - t107 * g(2);
t128 = -t106 * g(3) + t113 * t97;
t78 = (t106 * t113 * t115 - qJDD(2)) * pkin(2) + t128;
t62 = -t105 * t78 - t112 * t77;
t86 = qJD(4) + t88;
t35 = (-t81 * t86 - t54) * pkin(3) + t62;
t63 = -t105 * t77 + t112 * t78;
t124 = -t113 * g(3) - t106 * t97;
t83 = (-t106 ^ 2 * t115 - qJD(2) ^ 2) * pkin(2) + t124;
t125 = -t104 * t83 + t111 * t63;
t73 = -t89 * qJD(3) - t112 * qJDD(2) - t105 * t92;
t72 = qJDD(4) + t73;
t82 = -t104 * t98 + t111 * t89;
t36 = (-t81 * t82 + t72) * pkin(3) + t125;
t27 = t103 * t35 + t110 * t36;
t123 = -t104 * t63 - t111 * t83;
t43 = (-t82 ^ 2 - t86 ^ 2) * pkin(3) - t123;
t132 = t102 * t43 + t109 * t27;
t67 = t103 * t86 + t110 * t82;
t38 = -t67 * qJD(5) - t103 * t54 + t110 * t72;
t37 = qJDD(6) - t38;
t79 = qJD(5) - t81;
t55 = -t102 * t67 + t109 * t79;
t56 = t102 * t79 + t109 * t67;
t19 = (t55 * t56 - t37) * pkin(4) + t132;
t26 = t103 * t36 - t110 * t35;
t66 = -t103 * t82 + t110 * t86;
t39 = t66 * qJD(5) + t103 * t72 + t110 * t54;
t53 = -t82 * qJD(4) - t104 * t74 - t111 * t87;
t51 = qJDD(5) - t53;
t31 = t55 * qJD(6) + t102 * t51 + t109 * t39;
t65 = qJD(6) - t66;
t20 = (t55 * t65 + t31) * pkin(4) + t26;
t44 = -t101 * t56 - t108 * t65;
t23 = t44 * qJD(7) - t101 * t37 + t108 * t31;
t30 = -t56 * qJD(6) - t102 * t39 + t109 * t51;
t29 = qJDD(7) + t30;
t45 = -t101 * t65 + t108 * t56;
t32 = -t44 * mrSges(8,1) + t45 * mrSges(8,2);
t52 = qJD(7) + t55;
t33 = -t52 * mrSges(8,2) + t44 * mrSges(8,3);
t17 = m(8) * (-t101 * t19 - t108 * t20) - t23 * mrSges(8,3) + t29 * mrSges(8,1) - t45 * t32 + t52 * t33;
t22 = -t45 * qJD(7) - t101 * t31 - t108 * t37;
t34 = t52 * mrSges(8,1) - t45 * mrSges(8,3);
t18 = m(8) * (-t101 * t20 + t108 * t19) + t22 * mrSges(8,3) - t29 * mrSges(8,2) + t44 * t32 - t52 * t34;
t40 = -t55 * mrSges(7,1) + t56 * mrSges(7,2);
t47 = t65 * mrSges(7,1) - t56 * mrSges(7,3);
t15 = m(7) * t132 - t37 * mrSges(7,2) + t30 * mrSges(7,3) - t101 * t17 + t108 * t18 + t55 * t40 - t65 * t47;
t126 = -t102 * t27 + t109 * t43;
t121 = -t22 * mrSges(8,1) + t23 * mrSges(8,2) - t44 * t33 + m(8) * ((-t56 ^ 2 - t65 ^ 2) * pkin(4) + t126) + t45 * t34;
t46 = -t65 * mrSges(7,2) + t55 * mrSges(7,3);
t16 = m(7) * t126 + t37 * mrSges(7,1) - t31 * mrSges(7,3) - t56 * t40 + t65 * t46 + t121;
t57 = -t79 * mrSges(6,2) + t66 * mrSges(6,3);
t58 = t79 * mrSges(6,1) - t67 * mrSges(6,3);
t116 = m(6) * t43 - t38 * mrSges(6,1) + t39 * mrSges(6,2) + t102 * t15 + t109 * t16 - t66 * t57 + t67 * t58;
t64 = -t81 * mrSges(5,1) + t82 * mrSges(5,2);
t68 = -t86 * mrSges(5,2) + t81 * mrSges(5,3);
t10 = m(5) * t123 + t72 * mrSges(5,1) - t54 * mrSges(5,3) - t82 * t64 + t86 * t68 - t116;
t76 = -t88 * mrSges(4,1) + t89 * mrSges(4,2);
t85 = t98 * mrSges(4,1) - t89 * mrSges(4,3);
t48 = -t66 * mrSges(6,1) + t67 * mrSges(6,2);
t12 = m(6) * t27 - t51 * mrSges(6,2) + t38 * mrSges(6,3) - t102 * t16 + t109 * t15 + t66 * t48 - t79 * t58;
t117 = t30 * mrSges(7,1) - t31 * mrSges(7,2) + t101 * t18 + t108 * t17 + t55 * t46 - t56 * t47;
t14 = t51 * mrSges(6,1) - t39 * mrSges(6,3) - t67 * t48 + t79 * t57 + (-m(6) - m(7)) * t26 + t117;
t69 = t86 * mrSges(5,1) - t82 * mrSges(5,3);
t9 = m(5) * t125 - t72 * mrSges(5,2) + t53 * mrSges(5,3) - t103 * t14 + t110 * t12 + t81 * t64 - t86 * t69;
t7 = m(4) * t63 - t87 * mrSges(4,2) + t73 * mrSges(4,3) - t104 * t10 + t111 * t9 + t88 * t76 - t98 * t85;
t120 = m(5) * t62 - t53 * mrSges(5,1) + t54 * mrSges(5,2) + t103 * t12 + t110 * t14 - t81 * t68 + t82 * t69;
t84 = -t98 * mrSges(4,2) + t88 * mrSges(4,3);
t8 = m(4) * t62 + t87 * mrSges(4,1) - t74 * mrSges(4,3) - t89 * t76 + t98 * t84 + t120;
t94 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t131;
t95 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t130;
t134 = t93 * mrSges(3,1) - t92 * mrSges(3,2) - (t106 * t94 - t113 * t95) * qJD(1) + t105 * t7 + t112 * t8;
t91 = (-mrSges(3,1) * t113 + mrSges(3,2) * t106) * qJD(1);
t4 = m(3) * t128 - qJDD(2) * mrSges(3,2) + t93 * mrSges(3,3) - qJD(2) * t94 - t105 * t8 + t112 * t7 + t91 * t130;
t118 = m(4) * t83 - t73 * mrSges(4,1) + t74 * mrSges(4,2) - t111 * t10 - t104 * t9 - t88 * t84 + t89 * t85;
t6 = m(3) * t124 + qJDD(2) * mrSges(3,1) - t92 * mrSges(3,3) + qJD(2) * t95 - t91 * t131 + t118;
t133 = t106 * t4 + t113 * t6;
t2 = qJDD(1) * mrSges(2,1) - t115 * mrSges(2,2) + (m(2) + m(3)) * t96 + t134;
t1 = m(2) * t97 - t115 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t106 * t6 + t113 * t4;
t3 = [-m(1) * g(1) + t114 * t1 - t107 * t2, t1, t4, t7, t9, t12, t15, t18; -m(1) * g(2) + t107 * t1 + t114 * t2, t2, t6, t8, t10, t14, t16, t17; (-m(1) - m(2)) * g(3) + t133, -m(2) * g(3) + t133, -m(3) * t96 - t134, t118, t120, t116, m(7) * t26 - t117, t121;];
f_new  = t3;
