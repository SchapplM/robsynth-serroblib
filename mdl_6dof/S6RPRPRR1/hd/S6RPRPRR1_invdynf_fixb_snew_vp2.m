% Calculate vector of cutting forces with Newton-Euler
% S6RPRPRR1
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2,theta4]';
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
% Datum: 2019-05-05 18:16
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6RPRPRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 18:13:26
% EndTime: 2019-05-05 18:13:32
% DurationCPUTime: 3.56s
% Computational Cost: add. (48234->179), mult. (106180->240), div. (0->0), fcn. (73941->12), ass. (0->94)
t100 = cos(qJ(3));
t102 = qJD(1) ^ 2;
t101 = cos(qJ(1));
t97 = sin(qJ(1));
t117 = t97 * g(1) - t101 * g(2);
t73 = qJDD(1) * pkin(1) + t117;
t112 = -t101 * g(1) - t97 * g(2);
t75 = -t102 * pkin(1) + t112;
t91 = sin(pkin(10));
t93 = cos(pkin(10));
t115 = t93 * t73 - t91 * t75;
t109 = -qJDD(1) * pkin(2) - t115;
t96 = sin(qJ(3));
t122 = qJD(1) * t96;
t120 = qJD(1) * qJD(3);
t77 = t100 * qJDD(1) - t96 * t120;
t78 = qJD(3) * pkin(3) - qJ(4) * t122;
t88 = t100 ^ 2;
t106 = -t77 * pkin(3) + qJDD(4) + t78 * t122 + (-qJ(4) * t88 - pkin(7)) * t102 + t109;
t114 = t100 * t120;
t76 = t96 * qJDD(1) + t114;
t90 = sin(pkin(11));
t92 = cos(pkin(11));
t59 = -t90 * t76 + t92 * t77;
t68 = (t100 * t90 + t92 * t96) * qJD(1);
t63 = qJD(3) * pkin(4) - t68 * pkin(8);
t67 = (t100 * t92 - t90 * t96) * qJD(1);
t66 = t67 ^ 2;
t104 = -t59 * pkin(4) - t66 * pkin(8) + t68 * t63 + t106;
t123 = t91 * t73 + t93 * t75;
t54 = -t102 * pkin(2) + qJDD(1) * pkin(7) + t123;
t89 = -g(3) + qJDD(2);
t116 = t100 * t89 - t96 * t54;
t41 = (-t76 + t114) * qJ(4) + (t100 * t102 * t96 + qJDD(3)) * pkin(3) + t116;
t124 = t100 * t54 + t96 * t89;
t42 = -t88 * t102 * pkin(3) + t77 * qJ(4) - qJD(3) * t78 + t124;
t113 = -0.2e1 * qJD(4) * t68 + t92 * t41 - t90 * t42;
t60 = t92 * t76 + t90 * t77;
t21 = (qJD(3) * t67 - t60) * pkin(8) + (t67 * t68 + qJDD(3)) * pkin(4) + t113;
t118 = 0.2e1 * qJD(4) * t67 + t90 * t41 + t92 * t42;
t23 = -t66 * pkin(4) + t59 * pkin(8) - qJD(3) * t63 + t118;
t95 = sin(qJ(5));
t99 = cos(qJ(5));
t125 = t95 * t21 + t99 * t23;
t51 = t99 * t67 - t95 * t68;
t52 = t95 * t67 + t99 * t68;
t40 = -t51 * pkin(5) - t52 * pkin(9);
t87 = qJD(3) + qJD(5);
t85 = t87 ^ 2;
t86 = qJDD(3) + qJDD(5);
t18 = -t85 * pkin(5) + t86 * pkin(9) + t51 * t40 + t125;
t30 = -t52 * qJD(5) + t99 * t59 - t95 * t60;
t31 = t51 * qJD(5) + t95 * t59 + t99 * t60;
t19 = (-t51 * t87 - t31) * pkin(9) + (t52 * t87 - t30) * pkin(5) + t104;
t94 = sin(qJ(6));
t98 = cos(qJ(6));
t45 = -t94 * t52 + t98 * t87;
t25 = t45 * qJD(6) + t98 * t31 + t94 * t86;
t29 = qJDD(6) - t30;
t46 = t98 * t52 + t94 * t87;
t32 = -t45 * mrSges(7,1) + t46 * mrSges(7,2);
t50 = qJD(6) - t51;
t33 = -t50 * mrSges(7,2) + t45 * mrSges(7,3);
t15 = m(7) * (-t94 * t18 + t98 * t19) - t25 * mrSges(7,3) + t29 * mrSges(7,1) - t46 * t32 + t50 * t33;
t24 = -t46 * qJD(6) - t94 * t31 + t98 * t86;
t34 = t50 * mrSges(7,1) - t46 * mrSges(7,3);
t16 = m(7) * (t98 * t18 + t94 * t19) + t24 * mrSges(7,3) - t29 * mrSges(7,2) + t45 * t32 - t50 * t34;
t47 = -t87 * mrSges(6,2) + t51 * mrSges(6,3);
t48 = t87 * mrSges(6,1) - t52 * mrSges(6,3);
t108 = -m(6) * t104 + t30 * mrSges(6,1) - t31 * mrSges(6,2) - t98 * t15 - t94 * t16 + t51 * t47 - t52 * t48;
t61 = -qJD(3) * mrSges(5,2) + t67 * mrSges(5,3);
t62 = qJD(3) * mrSges(5,1) - t68 * mrSges(5,3);
t105 = -m(5) * t106 + t59 * mrSges(5,1) - t60 * mrSges(5,2) + t67 * t61 - t68 * t62 + t108;
t79 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t122;
t121 = qJD(1) * t100;
t80 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t121;
t126 = -(t100 * t80 - t96 * t79) * qJD(1) + m(4) * (-t102 * pkin(7) + t109) - t77 * mrSges(4,1) + t76 * mrSges(4,2) - t105;
t74 = (-mrSges(4,1) * t100 + mrSges(4,2) * t96) * qJD(1);
t39 = -t51 * mrSges(6,1) + t52 * mrSges(6,2);
t11 = m(6) * t125 - t86 * mrSges(6,2) + t30 * mrSges(6,3) - t94 * t15 + t98 * t16 + t51 * t39 - t87 * t48;
t111 = t99 * t21 - t95 * t23;
t107 = m(7) * (-t86 * pkin(5) - t85 * pkin(9) + t52 * t40 - t111) - t24 * mrSges(7,1) + t25 * mrSges(7,2) - t45 * t33 + t46 * t34;
t12 = m(6) * t111 + t86 * mrSges(6,1) - t31 * mrSges(6,3) - t52 * t39 + t87 * t47 - t107;
t57 = -t67 * mrSges(5,1) + t68 * mrSges(5,2);
t8 = m(5) * t113 + qJDD(3) * mrSges(5,1) - t60 * mrSges(5,3) + qJD(3) * t61 + t95 * t11 + t99 * t12 - t68 * t57;
t9 = m(5) * t118 - qJDD(3) * mrSges(5,2) + t59 * mrSges(5,3) - qJD(3) * t62 + t99 * t11 - t95 * t12 + t67 * t57;
t6 = m(4) * t116 + qJDD(3) * mrSges(4,1) - t76 * mrSges(4,3) + qJD(3) * t80 - t74 * t122 + t92 * t8 + t90 * t9;
t7 = m(4) * t124 - qJDD(3) * mrSges(4,2) + t77 * mrSges(4,3) - qJD(3) * t79 + t74 * t121 - t90 * t8 + t92 * t9;
t119 = m(3) * t89 + t100 * t6 + t96 * t7;
t10 = m(3) * t115 + qJDD(1) * mrSges(3,1) - t102 * mrSges(3,2) - t126;
t3 = m(3) * t123 - t102 * mrSges(3,1) - qJDD(1) * mrSges(3,2) + t100 * t7 - t96 * t6;
t2 = m(2) * t112 - t102 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t91 * t10 + t93 * t3;
t1 = m(2) * t117 + qJDD(1) * mrSges(2,1) - t102 * mrSges(2,2) + t93 * t10 + t91 * t3;
t4 = [-m(1) * g(1) - t97 * t1 + t101 * t2, t2, t3, t7, t9, t11, t16; -m(1) * g(2) + t101 * t1 + t97 * t2, t1, t10, t6, t8, t12, t15; (-m(1) - m(2)) * g(3) + t119, -m(2) * g(3) + t119, t119, t126, -t105, -t108, t107;];
f_new  = t4;
