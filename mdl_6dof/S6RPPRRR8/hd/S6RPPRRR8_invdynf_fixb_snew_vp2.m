% Calculate vector of cutting forces with Newton-Euler
% S6RPPRRR8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta3]';
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
% Datum: 2019-05-05 16:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPPRRR8_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR8_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 16:13:48
% EndTime: 2019-05-05 16:13:53
% DurationCPUTime: 2.27s
% Computational Cost: add. (28066->168), mult. (63657->212), div. (0->0), fcn. (45317->10), ass. (0->92)
t93 = sin(qJ(1));
t97 = cos(qJ(1));
t117 = t93 * g(1) - t97 * g(2);
t99 = qJD(1) ^ 2;
t107 = -t99 * qJ(2) + qJDD(2) - t117;
t128 = -pkin(1) - qJ(3);
t134 = -(2 * qJD(1) * qJD(3)) + t128 * qJDD(1) + t107;
t88 = sin(pkin(10));
t85 = t88 ^ 2;
t89 = cos(pkin(10));
t125 = t89 ^ 2 + t85;
t116 = t125 * mrSges(4,3);
t112 = -t97 * g(1) - t93 * g(2);
t133 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t112;
t132 = -m(2) - m(3);
t131 = pkin(3) * t99;
t130 = mrSges(2,1) - mrSges(3,2);
t129 = -mrSges(2,2) + mrSges(3,3);
t120 = t88 * g(3) + t134 * t89;
t40 = (-pkin(7) * qJDD(1) - t88 * t131) * t89 + t120;
t113 = -t89 * g(3) + t134 * t88;
t123 = qJDD(1) * t88;
t41 = -pkin(7) * t123 - t85 * t131 + t113;
t92 = sin(qJ(4));
t96 = cos(qJ(4));
t126 = t92 * t40 + t96 * t41;
t70 = (-t88 * t96 - t89 * t92) * qJD(1);
t111 = -t88 * t92 + t89 * t96;
t71 = t111 * qJD(1);
t52 = -t70 * pkin(4) - t71 * pkin(8);
t98 = qJD(4) ^ 2;
t23 = -t98 * pkin(4) + qJDD(4) * pkin(8) + t70 * t52 + t126;
t105 = qJDD(3) + t133;
t101 = pkin(3) * t123 + (-t125 * pkin(7) + t128) * t99 + t105;
t124 = t70 * qJD(4);
t122 = qJDD(1) * t89;
t67 = t71 * qJD(4);
t53 = -t92 * t122 - t96 * t123 - t67;
t54 = t111 * qJDD(1) + t124;
t26 = (-t54 - t124) * pkin(8) + (-t53 + t67) * pkin(4) + t101;
t91 = sin(qJ(5));
t95 = cos(qJ(5));
t127 = t95 * t23 + t91 * t26;
t68 = qJD(5) - t70;
t109 = -qJDD(1) * mrSges(4,3) - t99 * (mrSges(4,1) * t88 + mrSges(4,2) * t89);
t58 = t91 * qJD(4) + t95 * t71;
t31 = -t58 * qJD(5) + t95 * qJDD(4) - t91 * t54;
t57 = t95 * qJD(4) - t91 * t71;
t32 = t57 * qJD(5) + t91 * qJDD(4) + t95 * t54;
t90 = sin(qJ(6));
t94 = cos(qJ(6));
t34 = t90 * t57 + t94 * t58;
t19 = -t34 * qJD(6) + t94 * t31 - t90 * t32;
t33 = t94 * t57 - t90 * t58;
t20 = t33 * qJD(6) + t90 * t31 + t94 * t32;
t114 = t96 * t40 - t92 * t41;
t22 = -qJDD(4) * pkin(4) - t98 * pkin(8) + t71 * t52 - t114;
t66 = qJD(6) + t68;
t29 = -t66 * mrSges(7,2) + t33 * mrSges(7,3);
t30 = t66 * mrSges(7,1) - t34 * mrSges(7,3);
t44 = t68 * pkin(5) - t58 * pkin(9);
t56 = t57 ^ 2;
t106 = t19 * mrSges(7,1) + t33 * t29 - m(7) * (-t31 * pkin(5) - t56 * pkin(9) + t58 * t44 + t22) - t20 * mrSges(7,2) - t34 * t30;
t42 = -t68 * mrSges(6,2) + t57 * mrSges(6,3);
t43 = t68 * mrSges(6,1) - t58 * mrSges(6,3);
t100 = m(6) * t22 - t31 * mrSges(6,1) + t32 * mrSges(6,2) - t57 * t42 + t58 * t43 - t106;
t49 = -t70 * mrSges(5,1) + t71 * mrSges(5,2);
t62 = -qJD(4) * mrSges(5,2) + t70 * mrSges(5,3);
t11 = m(5) * t114 + qJDD(4) * mrSges(5,1) - t54 * mrSges(5,3) + qJD(4) * t62 - t71 * t49 - t100;
t115 = -t91 * t23 + t95 * t26;
t51 = qJDD(5) - t53;
t14 = (t57 * t68 - t32) * pkin(9) + (t57 * t58 + t51) * pkin(5) + t115;
t15 = -t56 * pkin(5) + t31 * pkin(9) - t68 * t44 + t127;
t28 = -t33 * mrSges(7,1) + t34 * mrSges(7,2);
t48 = qJDD(6) + t51;
t12 = m(7) * (t94 * t14 - t90 * t15) - t20 * mrSges(7,3) + t48 * mrSges(7,1) - t34 * t28 + t66 * t29;
t13 = m(7) * (t90 * t14 + t94 * t15) + t19 * mrSges(7,3) - t48 * mrSges(7,2) + t33 * t28 - t66 * t30;
t36 = -t57 * mrSges(6,1) + t58 * mrSges(6,2);
t10 = m(6) * t127 - t51 * mrSges(6,2) + t31 * mrSges(6,3) - t90 * t12 + t94 * t13 + t57 * t36 - t68 * t43;
t63 = qJD(4) * mrSges(5,1) - t71 * mrSges(5,3);
t9 = m(6) * t115 + t51 * mrSges(6,1) - t32 * mrSges(6,3) + t94 * t12 + t90 * t13 - t58 * t36 + t68 * t42;
t6 = m(5) * t126 - qJDD(4) * mrSges(5,2) + t53 * mrSges(5,3) - qJD(4) * t63 + t95 * t10 + t70 * t49 - t91 * t9;
t3 = m(4) * t120 + t109 * t89 + t96 * t11 + t92 * t6;
t4 = m(4) * t113 + t109 * t88 - t92 * t11 + t96 * t6;
t118 = -t88 * t3 + t89 * t4;
t108 = -m(3) * (-qJDD(1) * pkin(1) + t107) - t89 * t3 - t88 * t4;
t104 = -m(5) * t101 + t53 * mrSges(5,1) - t54 * mrSges(5,2) - t91 * t10 + t70 * t62 - t71 * t63 - t95 * t9;
t103 = -t104 + m(4) * (t128 * t99 + t105) + mrSges(4,1) * t123 + mrSges(4,2) * t122;
t102 = m(3) * (t99 * pkin(1) - t133) - t103;
t5 = (-t116 - t130) * t99 + t129 * qJDD(1) + m(2) * t112 - t102;
t1 = m(2) * t117 + t130 * qJDD(1) + t129 * t99 + t108;
t2 = [-m(1) * g(1) - t93 * t1 + t97 * t5, t5, -m(3) * g(3) + t118, t4, t6, t10, t13; -m(1) * g(2) + t97 * t1 + t93 * t5, t1, t102 - qJDD(1) * mrSges(3,3) + (-mrSges(3,2) + t116) * t99, t3, t11, t9, t12; (-m(1) + t132) * g(3) + t118, t132 * g(3) + t118, qJDD(1) * mrSges(3,2) - t99 * mrSges(3,3) - t108, -t99 * t116 + t103, -t104, t100, -t106;];
f_new  = t2;
