% Calculate vector of cutting forces with Newton-Euler
% S6PPRRRR1
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
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d3,d4,d5,d6,theta1,theta2]';
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
% Datum: 2019-05-04 20:42
% Revision: 89c353f7eff3bd693eda4e29f35b2761dbc3ada0 (2019-05-03)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S6PPRRRR1_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(3,1),zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: qJD has to be [6x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [6 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: qJDD has to be [6x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PPRRRR1_invdynf_fixb_snew_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-05-04 20:39:04
% EndTime: 2019-05-04 20:39:10
% DurationCPUTime: 4.42s
% Computational Cost: add. (68012->147), mult. (125357->207), div. (0->0), fcn. (98638->16), ass. (0->91)
t77 = sin(pkin(12));
t81 = cos(pkin(12));
t65 = t77 * g(1) - t81 * g(2);
t75 = -g(3) + qJDD(1);
t79 = sin(pkin(6));
t83 = cos(pkin(6));
t100 = t65 * t83 + t75 * t79;
t66 = -t81 * g(1) - t77 * g(2);
t76 = sin(pkin(13));
t80 = cos(pkin(13));
t47 = t100 * t80 - t76 * t66;
t57 = -t79 * t65 + t83 * t75 + qJDD(2);
t78 = sin(pkin(7));
t82 = cos(pkin(7));
t120 = t47 * t82 + t57 * t78;
t48 = t100 * t76 + t80 * t66;
t87 = sin(qJ(3));
t91 = cos(qJ(3));
t107 = t120 * t87 + t91 * t48;
t92 = qJD(3) ^ 2;
t31 = -t92 * pkin(3) + qJDD(3) * pkin(9) + t107;
t36 = -t78 * t47 + t82 * t57;
t86 = sin(qJ(4));
t90 = cos(qJ(4));
t104 = -t86 * t31 + t90 * t36;
t110 = qJD(3) * t86;
t108 = qJD(3) * qJD(4);
t105 = t90 * t108;
t63 = t86 * qJDD(3) + t105;
t25 = (-t63 + t105) * pkin(10) + (t86 * t90 * t92 + qJDD(4)) * pkin(4) + t104;
t111 = t90 * t31 + t86 * t36;
t64 = t90 * qJDD(3) - t86 * t108;
t69 = qJD(4) * pkin(4) - pkin(10) * t110;
t74 = t90 ^ 2;
t26 = -t74 * t92 * pkin(4) + t64 * pkin(10) - qJD(4) * t69 + t111;
t85 = sin(qJ(5));
t89 = cos(qJ(5));
t112 = t85 * t25 + t89 * t26;
t59 = (t85 * t86 - t89 * t90) * qJD(3);
t60 = (t85 * t90 + t86 * t89) * qJD(3);
t51 = t59 * pkin(5) - t60 * pkin(11);
t73 = qJD(4) + qJD(5);
t71 = t73 ^ 2;
t72 = qJDD(4) + qJDD(5);
t21 = -t71 * pkin(5) + t72 * pkin(11) - t59 * t51 + t112;
t40 = -t60 * qJD(5) - t85 * t63 + t89 * t64;
t41 = -t59 * qJD(5) + t89 * t63 + t85 * t64;
t118 = t120 * t91 - t87 * t48;
t96 = -qJDD(3) * pkin(3) - t118;
t93 = -t64 * pkin(4) + t69 * t110 + (-pkin(10) * t74 - pkin(9)) * t92 + t96;
t22 = (t59 * t73 - t41) * pkin(11) + (t60 * t73 - t40) * pkin(5) + t93;
t84 = sin(qJ(6));
t88 = cos(qJ(6));
t53 = -t84 * t60 + t88 * t73;
t33 = t53 * qJD(6) + t88 * t41 + t84 * t72;
t54 = t88 * t60 + t84 * t73;
t37 = -t53 * mrSges(7,1) + t54 * mrSges(7,2);
t39 = qJDD(6) - t40;
t58 = qJD(6) + t59;
t42 = -t58 * mrSges(7,2) + t53 * mrSges(7,3);
t18 = m(7) * (-t84 * t21 + t88 * t22) - t33 * mrSges(7,3) + t39 * mrSges(7,1) - t54 * t37 + t58 * t42;
t32 = -t54 * qJD(6) - t84 * t41 + t88 * t72;
t43 = t58 * mrSges(7,1) - t54 * mrSges(7,3);
t19 = m(7) * (t88 * t21 + t84 * t22) + t32 * mrSges(7,3) - t39 * mrSges(7,2) + t53 * t37 - t58 * t43;
t50 = t59 * mrSges(6,1) + t60 * mrSges(6,2);
t56 = t73 * mrSges(6,1) - t60 * mrSges(6,3);
t14 = m(6) * t112 - t72 * mrSges(6,2) + t40 * mrSges(6,3) - t84 * t18 + t88 * t19 - t59 * t50 - t73 * t56;
t102 = t89 * t25 - t85 * t26;
t55 = -t73 * mrSges(6,2) - t59 * mrSges(6,3);
t95 = m(7) * (-t72 * pkin(5) - t71 * pkin(11) + t60 * t51 - t102) - t32 * mrSges(7,1) + t33 * mrSges(7,2) - t53 * t42 + t54 * t43;
t15 = m(6) * t102 + t72 * mrSges(6,1) - t41 * mrSges(6,3) - t60 * t50 + t73 * t55 - t95;
t62 = (-mrSges(5,1) * t90 + mrSges(5,2) * t86) * qJD(3);
t109 = qJD(3) * t90;
t68 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t109;
t11 = m(5) * t104 + qJDD(4) * mrSges(5,1) - t63 * mrSges(5,3) + qJD(4) * t68 - t62 * t110 + t85 * t14 + t89 * t15;
t67 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t110;
t12 = m(5) * t111 - qJDD(4) * mrSges(5,2) + t64 * mrSges(5,3) - qJD(4) * t67 + t62 * t109 + t89 * t14 - t85 * t15;
t10 = m(4) * t36 + t90 * t11 + t86 * t12;
t97 = -m(6) * t93 + t40 * mrSges(6,1) - t41 * mrSges(6,2) - t88 * t18 - t84 * t19 - t59 * t55 - t60 * t56;
t117 = (t86 * t67 - t90 * t68) * qJD(3) + m(5) * (-t92 * pkin(9) + t96) - t64 * mrSges(5,1) + t63 * mrSges(5,2) - t97;
t13 = m(4) * t118 + qJDD(3) * mrSges(4,1) - t92 * mrSges(4,2) - t117;
t9 = m(4) * t107 - t92 * mrSges(4,1) - qJDD(3) * mrSges(4,2) - t86 * t11 + t90 * t12;
t103 = t13 * t91 + t87 * t9;
t4 = m(3) * t47 - t78 * t10 + t103 * t82;
t8 = m(3) * t48 - t87 * t13 + t91 * t9;
t119 = t4 * t80 + t76 * t8;
t6 = m(3) * t57 + t82 * t10 + t103 * t78;
t106 = m(2) * t75 + t119 * t79 + t83 * t6;
t2 = m(2) * t66 - t76 * t4 + t80 * t8;
t1 = m(2) * t65 + t119 * t83 - t79 * t6;
t3 = [-m(1) * g(1) - t77 * t1 + t81 * t2, t2, t8, t9, t12, t14, t19; -m(1) * g(2) + t81 * t1 + t77 * t2, t1, t4, t13, t11, t15, t18; -m(1) * g(3) + t106, t106, t6, t10, t117, -t97, t95;];
f_new  = t3;
