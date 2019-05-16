% Calculate vector of cutting forces with Newton-Euler
% S6PRPRRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,d6,theta1]';
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
% Datum: 2019-05-05 01:35
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PRPRRR6_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRR6_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-05 01:30:43
% EndTime: 2019-05-05 01:30:48
% DurationCPUTime: 1.73s
% Computational Cost: add. (22874->150), mult. (42901->196), div. (0->0), fcn. (28542->12), ass. (0->89)
t84 = sin(pkin(6));
t90 = sin(qJ(2));
t119 = t84 * t90;
t83 = sin(pkin(11));
t85 = cos(pkin(11));
t69 = t83 * g(1) - t85 * g(2);
t86 = cos(pkin(6));
t120 = t69 * t86;
t70 = -t85 * g(1) - t83 * g(2);
t82 = -g(3) + qJDD(1);
t94 = cos(qJ(2));
t111 = t82 * t119 + t90 * t120 + t94 * t70;
t124 = -qJDD(2) * qJ(3) - (2 * qJD(3) * qJD(2)) - t111;
t123 = (t82 * t84 + t120) * t94 - t90 * t70;
t122 = -pkin(2) - pkin(8);
t96 = qJD(2) ^ 2;
t98 = -t96 * qJ(3) + qJDD(3) - t123;
t36 = t122 * qJDD(2) + t98;
t52 = -t84 * t69 + t86 * t82;
t89 = sin(qJ(4));
t93 = cos(qJ(4));
t115 = t89 * t36 + t93 * t52;
t66 = (pkin(4) * t89 - pkin(9) * t93) * qJD(2);
t78 = t89 * qJD(2);
t95 = qJD(4) ^ 2;
t27 = -t95 * pkin(4) + qJDD(4) * pkin(9) - t66 * t78 + t115;
t100 = t122 * t96 - t124;
t113 = qJD(2) * qJD(4);
t109 = t89 * t113;
t76 = t93 * t113;
t67 = -t89 * qJDD(2) - t76;
t68 = t93 * qJDD(2) - t109;
t30 = (-t68 + t109) * pkin(9) + (-t67 + t76) * pkin(4) + t100;
t88 = sin(qJ(5));
t92 = cos(qJ(5));
t108 = -t88 * t27 + t92 * t30;
t114 = qJD(2) * t93;
t63 = t92 * qJD(4) - t88 * t114;
t42 = t63 * qJD(5) + t88 * qJDD(4) + t92 * t68;
t60 = qJDD(5) - t67;
t64 = t88 * qJD(4) + t92 * t114;
t75 = t78 + qJD(5);
t18 = (t63 * t75 - t42) * pkin(10) + (t63 * t64 + t60) * pkin(5) + t108;
t116 = t92 * t27 + t88 * t30;
t41 = -t64 * qJD(5) + t92 * qJDD(4) - t88 * t68;
t50 = t75 * pkin(5) - t64 * pkin(10);
t59 = t63 ^ 2;
t19 = -t59 * pkin(5) + t41 * pkin(10) - t75 * t50 + t116;
t87 = sin(qJ(6));
t91 = cos(qJ(6));
t43 = t91 * t63 - t87 * t64;
t24 = t43 * qJD(6) + t87 * t41 + t91 * t42;
t44 = t87 * t63 + t91 * t64;
t32 = -t43 * mrSges(7,1) + t44 * mrSges(7,2);
t74 = qJD(6) + t75;
t39 = -t74 * mrSges(7,2) + t43 * mrSges(7,3);
t55 = qJDD(6) + t60;
t16 = m(7) * (t91 * t18 - t87 * t19) - t24 * mrSges(7,3) + t55 * mrSges(7,1) - t44 * t32 + t74 * t39;
t23 = -t44 * qJD(6) + t91 * t41 - t87 * t42;
t40 = t74 * mrSges(7,1) - t44 * mrSges(7,3);
t17 = m(7) * (t87 * t18 + t91 * t19) + t23 * mrSges(7,3) - t55 * mrSges(7,2) + t43 * t32 - t74 * t40;
t45 = -t63 * mrSges(6,1) + t64 * mrSges(6,2);
t48 = -t75 * mrSges(6,2) + t63 * mrSges(6,3);
t13 = m(6) * t108 + t60 * mrSges(6,1) - t42 * mrSges(6,3) + t91 * t16 + t87 * t17 - t64 * t45 + t75 * t48;
t49 = t75 * mrSges(6,1) - t64 * mrSges(6,3);
t14 = m(6) * t116 - t60 * mrSges(6,2) + t41 * mrSges(6,3) - t87 * t16 + t91 * t17 + t63 * t45 - t75 * t49;
t65 = (mrSges(5,1) * t89 + mrSges(5,2) * t93) * qJD(2);
t72 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t114;
t10 = m(5) * t115 - qJDD(4) * mrSges(5,2) + t67 * mrSges(5,3) - qJD(4) * t72 - t88 * t13 + t92 * t14 - t65 * t78;
t107 = t93 * t36 - t89 * t52;
t71 = -qJD(4) * mrSges(5,2) - mrSges(5,3) * t78;
t26 = -qJDD(4) * pkin(4) - t95 * pkin(9) + t66 * t114 - t107;
t102 = t23 * mrSges(7,1) + t43 * t39 - m(7) * (-t41 * pkin(5) - t59 * pkin(10) + t64 * t50 + t26) - t24 * mrSges(7,2) - t44 * t40;
t97 = m(6) * t26 - t41 * mrSges(6,1) + t42 * mrSges(6,2) - t63 * t48 + t64 * t49 - t102;
t15 = m(5) * t107 + qJDD(4) * mrSges(5,1) - t68 * mrSges(5,3) + qJD(4) * t71 - t65 * t114 - t97;
t103 = -m(4) * (-qJDD(2) * pkin(2) + t98) - t89 * t10 - t93 * t15;
t117 = (-mrSges(3,2) + mrSges(4,3));
t118 = mrSges(3,1) - mrSges(4,2);
t4 = m(3) * t123 + t118 * qJDD(2) + (t117 * t96) + t103;
t121 = t4 * t94;
t105 = m(4) * t52 + t93 * t10 - t89 * t15;
t6 = m(3) * t52 + t105;
t101 = m(5) * t100 - t67 * mrSges(5,1) + t68 * mrSges(5,2) + t72 * t114 + t92 * t13 + t88 * t14 + t71 * t78;
t99 = -m(4) * (t96 * pkin(2) + t124) + t101;
t8 = m(3) * t111 + t117 * qJDD(2) - t118 * t96 + t99;
t110 = m(2) * t82 + t8 * t119 + t84 * t121 + t86 * t6;
t2 = m(2) * t70 - t90 * t4 + t94 * t8;
t1 = m(2) * t69 - t84 * t6 + (t8 * t90 + t121) * t86;
t3 = [-m(1) * g(1) - t83 * t1 + t85 * t2, t2, t8, t105, t10, t14, t17; -m(1) * g(2) + t85 * t1 + t83 * t2, t1, t4, -(t96 * mrSges(4,2)) - qJDD(2) * mrSges(4,3) - t99, t15, t13, t16; -m(1) * g(3) + t110, t110, t6, qJDD(2) * mrSges(4,2) - t96 * mrSges(4,3) - t103, t101, t97, -t102;];
f_new  = t3;
