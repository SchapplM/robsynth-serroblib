% Calculate vector of cutting forces with Newton-Euler
% S6RPPRPR7
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
% Datum: 2019-05-05 14:33
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 14:31:14
% EndTime: 2019-05-05 14:31:18
% DurationCPUTime: 2.18s
% Computational Cost: add. (26199->165), mult. (61612->212), div. (0->0), fcn. (43355->10), ass. (0->91)
t89 = sin(qJ(1));
t92 = cos(qJ(1));
t113 = t89 * g(1) - t92 * g(2);
t94 = qJD(1) ^ 2;
t102 = -t94 * qJ(2) + qJDD(2) - t113;
t124 = -pkin(1) - qJ(3);
t131 = -(2 * qJD(1) * qJD(3)) + t124 * qJDD(1) + t102;
t84 = sin(pkin(9));
t80 = t84 ^ 2;
t86 = cos(pkin(9));
t122 = t86 ^ 2 + t80;
t112 = t122 * mrSges(4,3);
t108 = -t92 * g(1) - t89 * g(2);
t130 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t108;
t129 = -m(2) - m(3);
t128 = pkin(3) * t94;
t127 = mrSges(4,2) * t86;
t126 = mrSges(2,1) - mrSges(3,2);
t125 = -mrSges(2,2) + mrSges(3,3);
t116 = t84 * g(3) + t131 * t86;
t38 = (-pkin(7) * qJDD(1) - t84 * t128) * t86 + t116;
t110 = -t86 * g(3) + t131 * t84;
t119 = qJDD(1) * t84;
t39 = -pkin(7) * t119 - t80 * t128 + t110;
t88 = sin(qJ(4));
t91 = cos(qJ(4));
t123 = t88 * t38 + t91 * t39;
t107 = t84 * t91 + t86 * t88;
t69 = t107 * qJD(1);
t121 = t69 * qJD(4);
t106 = -t84 * t88 + t86 * t91;
t70 = t106 * qJD(1);
t120 = t70 * qJD(4);
t48 = t69 * pkin(4) - t70 * qJ(5);
t93 = qJD(4) ^ 2;
t20 = -t93 * pkin(4) + qJDD(4) * qJ(5) - t69 * t48 + t123;
t52 = t107 * qJDD(1) + t120;
t53 = t106 * qJDD(1) - t121;
t100 = qJDD(3) + t130;
t96 = pkin(3) * t119 + (-t122 * pkin(7) + t124) * t94 + t100;
t26 = (-t53 + t121) * qJ(5) + (t52 + t120) * pkin(4) + t96;
t83 = sin(pkin(10));
t85 = cos(pkin(10));
t57 = t85 * qJD(4) - t83 * t70;
t117 = 0.2e1 * qJD(5) * t57 + t85 * t20 + t83 * t26;
t104 = -qJDD(1) * mrSges(4,3) - t94 * (mrSges(4,1) * t84 + t127);
t111 = t91 * t38 - t88 * t39;
t49 = t69 * mrSges(5,1) + t70 * mrSges(5,2);
t63 = -qJD(4) * mrSges(5,2) - t69 * mrSges(5,3);
t19 = -qJDD(4) * pkin(4) - t93 * qJ(5) + t70 * t48 + qJDD(5) - t111;
t58 = t83 * qJD(4) + t85 * t70;
t87 = sin(qJ(6));
t90 = cos(qJ(6));
t32 = t87 * t57 + t90 * t58;
t43 = t85 * qJDD(4) - t83 * t53;
t44 = t83 * qJDD(4) + t85 * t53;
t24 = -t32 * qJD(6) + t90 * t43 - t87 * t44;
t31 = t90 * t57 - t87 * t58;
t25 = t31 * qJD(6) + t87 * t43 + t90 * t44;
t67 = qJD(6) + t69;
t29 = -t67 * mrSges(7,2) + t31 * mrSges(7,3);
t30 = t67 * mrSges(7,1) - t32 * mrSges(7,3);
t42 = t69 * pkin(5) - t58 * pkin(8);
t56 = t57 ^ 2;
t101 = t24 * mrSges(7,1) + t31 * t29 - m(7) * (-t43 * pkin(5) - t56 * pkin(8) + t58 * t42 + t19) - t25 * mrSges(7,2) - t32 * t30;
t40 = -t69 * mrSges(6,2) + t57 * mrSges(6,3);
t41 = t69 * mrSges(6,1) - t58 * mrSges(6,3);
t95 = m(6) * t19 - t43 * mrSges(6,1) + t44 * mrSges(6,2) - t57 * t40 + t58 * t41 - t101;
t11 = m(5) * t111 + qJDD(4) * mrSges(5,1) - t53 * mrSges(5,3) + qJD(4) * t63 - t70 * t49 - t95;
t109 = -0.2e1 * qJD(5) * t58 - t83 * t20 + t85 * t26;
t14 = (t57 * t69 - t44) * pkin(8) + (t57 * t58 + t52) * pkin(5) + t109;
t15 = -t56 * pkin(5) + t43 * pkin(8) - t69 * t42 + t117;
t28 = -t31 * mrSges(7,1) + t32 * mrSges(7,2);
t51 = qJDD(6) + t52;
t12 = m(7) * (t90 * t14 - t87 * t15) - t25 * mrSges(7,3) + t51 * mrSges(7,1) - t32 * t28 + t67 * t29;
t13 = m(7) * (t87 * t14 + t90 * t15) + t24 * mrSges(7,3) - t51 * mrSges(7,2) + t31 * t28 - t67 * t30;
t33 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t10 = m(6) * t117 - t52 * mrSges(6,2) + t43 * mrSges(6,3) - t87 * t12 + t90 * t13 + t57 * t33 - t69 * t41;
t64 = qJD(4) * mrSges(5,1) - t70 * mrSges(5,3);
t9 = m(6) * t109 + t52 * mrSges(6,1) - t44 * mrSges(6,3) + t90 * t12 + t87 * t13 - t58 * t33 + t69 * t40;
t6 = m(5) * t123 - qJDD(4) * mrSges(5,2) - t52 * mrSges(5,3) - qJD(4) * t64 + t85 * t10 - t69 * t49 - t83 * t9;
t3 = m(4) * t116 + t104 * t86 + t91 * t11 + t88 * t6;
t4 = m(4) * t110 + t104 * t84 - t88 * t11 + t91 * t6;
t114 = -t84 * t3 + t86 * t4;
t103 = -m(3) * (-qJDD(1) * pkin(1) + t102) - t86 * t3 - t84 * t4;
t99 = m(5) * t96 + t52 * mrSges(5,1) + t53 * mrSges(5,2) + t83 * t10 + t69 * t63 + t70 * t64 + t85 * t9;
t98 = m(4) * (t124 * t94 + t100) + mrSges(4,1) * t119 + qJDD(1) * t127 + t99;
t97 = m(3) * (t94 * pkin(1) - t130) - t98;
t5 = (-t112 - t126) * t94 + t125 * qJDD(1) + m(2) * t108 - t97;
t1 = m(2) * t113 + t126 * qJDD(1) + t125 * t94 + t103;
t2 = [-m(1) * g(1) - t89 * t1 + t92 * t5, t5, -m(3) * g(3) + t114, t4, t6, t10, t13; -m(1) * g(2) + t92 * t1 + t89 * t5, t1, t97 + (-mrSges(3,2) + t112) * t94 - qJDD(1) * mrSges(3,3), t3, t11, t9, t12; (-m(1) + t129) * g(3) + t114, t129 * g(3) + t114, qJDD(1) * mrSges(3,2) - t94 * mrSges(3,3) - t103, -t94 * t112 + t98, t99, t95, -t101;];
f_new  = t2;
