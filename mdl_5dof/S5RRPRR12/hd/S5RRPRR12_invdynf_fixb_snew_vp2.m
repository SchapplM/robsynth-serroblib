% Calculate vector of cutting forces with Newton-Euler
% S5RRPRR12
% Use Code from Maple symbolic Code Generation
%
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% qJDD [5x1]
%   Generalized joint accelerations
% g [3x1]
%   gravitation vector in mdh base frame [m/s^2]
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
%
% Output:
% f_new [3x6]
%   vector of cutting forces (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 20:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RRPRR12_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR12_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 20:29:04
% EndTime: 2019-12-31 20:29:07
% DurationCPUTime: 1.02s
% Computational Cost: add. (8460->162), mult. (17392->208), div. (0->0), fcn. (10237->8), ass. (0->80)
t91 = cos(qJ(2));
t111 = qJD(1) * t91;
t68 = mrSges(4,2) * t111 + qJD(2) * mrSges(4,3);
t114 = -qJD(2) * mrSges(3,2) + mrSges(3,3) * t111 + t68;
t88 = sin(qJ(1));
t92 = cos(qJ(1));
t113 = t88 * g(1) - t92 * g(2);
t94 = qJD(1) ^ 2;
t53 = -qJDD(1) * pkin(1) - t94 * pkin(6) - t113;
t110 = qJD(1) * qJD(2);
t107 = t91 * t110;
t87 = sin(qJ(2));
t62 = t87 * qJDD(1) + t107;
t108 = t87 * t110;
t63 = t91 * qJDD(1) - t108;
t112 = qJD(1) * t87;
t65 = qJD(2) * mrSges(3,1) - mrSges(3,3) * t112;
t66 = -qJD(2) * mrSges(4,1) + mrSges(4,2) * t112;
t102 = -t63 * pkin(2) + t53 + (-t107 - t62) * qJ(3);
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t52 = (-t86 * t91 + t87 * t90) * qJD(1);
t33 = -t52 * qJD(4) - t86 * t62 - t90 * t63;
t51 = (t86 * t87 + t90 * t91) * qJD(1);
t34 = -t51 * qJD(4) + t90 * t62 - t86 * t63;
t79 = -qJD(2) + qJD(4);
t120 = t91 ^ 2 * t94;
t122 = 2 * qJD(3);
t69 = -qJD(2) * pkin(3) - pkin(7) * t112;
t95 = -pkin(2) * t108 + t63 * pkin(3) - pkin(7) * t120 - t102 + (t122 + t69) * t112;
t14 = t95 + (t51 * t79 - t34) * pkin(8) + (t52 * t79 - t33) * pkin(4);
t106 = -t92 * g(1) - t88 * g(2);
t54 = -t94 * pkin(1) + qJDD(1) * pkin(6) + t106;
t109 = -t87 * g(3) + t91 * t54;
t59 = (-pkin(2) * t91 - qJ(3) * t87) * qJD(1);
t93 = qJD(2) ^ 2;
t98 = -t93 * pkin(2) + qJDD(2) * qJ(3) + qJD(2) * t122 + t59 * t111 + t109;
t23 = -pkin(3) * t120 - t63 * pkin(7) + qJD(2) * t69 + t98;
t116 = -t91 * g(3) - t87 * t54;
t32 = -qJDD(2) * pkin(2) - t93 * qJ(3) + t59 * t112 + qJDD(3) - t116;
t24 = (-t62 + t107) * pkin(7) + (-t87 * t91 * t94 - qJDD(2)) * pkin(3) + t32;
t117 = t90 * t23 + t86 * t24;
t40 = t51 * pkin(4) - t52 * pkin(8);
t77 = t79 ^ 2;
t78 = -qJDD(2) + qJDD(4);
t16 = -t77 * pkin(4) + t78 * pkin(8) - t51 * t40 + t117;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t41 = -t85 * t52 + t89 * t79;
t20 = t41 * qJD(5) + t89 * t34 + t85 * t78;
t42 = t89 * t52 + t85 * t79;
t27 = -t41 * mrSges(6,1) + t42 * mrSges(6,2);
t31 = qJDD(5) - t33;
t50 = qJD(5) + t51;
t35 = -t50 * mrSges(6,2) + t41 * mrSges(6,3);
t12 = m(6) * (t89 * t14 - t85 * t16) - t20 * mrSges(6,3) + t31 * mrSges(6,1) - t42 * t27 + t50 * t35;
t19 = -t42 * qJD(5) - t85 * t34 + t89 * t78;
t36 = t50 * mrSges(6,1) - t42 * mrSges(6,3);
t13 = m(6) * (t85 * t14 + t89 * t16) + t19 * mrSges(6,3) - t31 * mrSges(6,2) + t41 * t27 - t50 * t36;
t43 = -t79 * mrSges(5,2) - t51 * mrSges(5,3);
t44 = t79 * mrSges(5,1) - t52 * mrSges(5,3);
t103 = m(5) * t95 - t33 * mrSges(5,1) + t34 * mrSges(5,2) + t89 * t12 + t85 * t13 + t51 * t43 + t52 * t44;
t99 = m(4) * ((pkin(2) * qJD(2) - (2 * qJD(3))) * t112 + t102) - t63 * mrSges(4,1) - t103;
t125 = -(t114 * t91 - (t65 - t66) * t87) * qJD(1) + m(3) * t53 - t63 * mrSges(3,1) + (mrSges(3,2) - mrSges(4,3)) * t62 + t99;
t60 = (-mrSges(4,1) * t91 - mrSges(4,3) * t87) * qJD(1);
t39 = t51 * mrSges(5,1) + t52 * mrSges(5,2);
t8 = m(5) * t117 - t78 * mrSges(5,2) + t33 * mrSges(5,3) - t85 * t12 + t89 * t13 - t51 * t39 - t79 * t44;
t105 = -t86 * t23 + t90 * t24;
t97 = m(6) * (-t78 * pkin(4) - t77 * pkin(8) + t52 * t40 - t105) - t19 * mrSges(6,1) + t20 * mrSges(6,2) - t41 * t35 + t42 * t36;
t9 = m(5) * t105 + t78 * mrSges(5,1) - t34 * mrSges(5,3) - t52 * t39 + t79 * t43 - t97;
t101 = m(4) * t98 + qJDD(2) * mrSges(4,3) + qJD(2) * t66 + t60 * t111 + t90 * t8 - t86 * t9;
t118 = mrSges(3,3) + mrSges(4,2);
t61 = (-mrSges(3,1) * t91 + mrSges(3,2) * t87) * qJD(1);
t4 = m(3) * t109 - qJDD(2) * mrSges(3,2) - qJD(2) * t65 + t61 * t111 + t118 * t63 + t101;
t100 = -m(4) * t32 - t86 * t8 - t90 * t9;
t5 = m(3) * t116 - t118 * t62 + (mrSges(3,1) + mrSges(4,1)) * qJDD(2) + t114 * qJD(2) + (-t60 - t61) * t112 + t100;
t121 = t87 * t4 + t91 * t5;
t6 = m(2) * t113 + qJDD(1) * mrSges(2,1) - t94 * mrSges(2,2) - t125;
t1 = m(2) * t106 - t94 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t91 * t4 - t87 * t5;
t2 = [-m(1) * g(1) + t92 * t1 - t88 * t6, t1, t4, t63 * mrSges(4,2) + t101, t8, t13; -m(1) * g(2) + t88 * t1 + t92 * t6, t6, t5, -t62 * mrSges(4,3) + (-t87 * t66 - t91 * t68) * qJD(1) + t99, t9, t12; (-m(1) - m(2)) * g(3) + t121, -m(2) * g(3) + t121, t125, -qJDD(2) * mrSges(4,1) + t62 * mrSges(4,2) - qJD(2) * t68 + t60 * t112 - t100, t103, t97;];
f_new = t2;
