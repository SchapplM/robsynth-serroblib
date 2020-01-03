% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR11
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:28
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR11_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR11_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:27:13
% EndTime: 2019-12-31 18:27:15
% DurationCPUTime: 1.09s
% Computational Cost: add. (8553->151), mult. (20646->186), div. (0->0), fcn. (13918->8), ass. (0->78)
t113 = cos(qJ(3));
t75 = sin(pkin(8));
t76 = cos(pkin(8));
t78 = sin(qJ(3));
t119 = -t76 * t113 + t75 * t78;
t58 = t119 * qJD(1);
t105 = t58 * qJD(3);
t90 = t113 * t75 + t76 * t78;
t48 = t90 * qJDD(1) - t105;
t79 = sin(qJ(1));
t81 = cos(qJ(1));
t108 = t79 * g(1) - t81 * g(2);
t100 = -qJDD(2) + t108;
t71 = t76 ^ 2;
t107 = t75 ^ 2 + t71;
t83 = qJD(1) ^ 2;
t86 = -(t107 * pkin(6) + qJ(2)) * t83 - (pkin(2) * t76 + pkin(1)) * qJDD(1) - t100;
t122 = (-t48 + t105) * qJ(4) + t86;
t121 = t107 * mrSges(3,3);
t116 = 2 * qJD(4);
t59 = t90 * qJD(1);
t42 = t58 * mrSges(5,1) - t59 * mrSges(5,3);
t109 = -t58 * mrSges(4,1) - t59 * mrSges(4,2) - t42;
t106 = pkin(6) * qJDD(1);
t114 = pkin(2) * t83;
t95 = -t81 * g(1) - t79 * g(2);
t60 = -t83 * pkin(1) + qJDD(1) * qJ(2) + t95;
t103 = qJD(1) * qJD(2);
t99 = -t76 * g(3) - 0.2e1 * t75 * t103;
t34 = (t76 * t114 - t106 - t60) * t75 + t99;
t97 = -t75 * g(3) + (0.2e1 * t103 + t60) * t76;
t35 = t76 * t106 - t71 * t114 + t97;
t110 = t113 * t35 + t78 * t34;
t111 = -mrSges(4,3) - mrSges(5,2);
t104 = t59 * qJD(3);
t47 = t119 * qJDD(1) + t104;
t51 = qJD(3) * mrSges(4,1) - t59 * mrSges(4,3);
t41 = t58 * pkin(3) - t59 * qJ(4);
t82 = qJD(3) ^ 2;
t94 = t113 * t34 - t78 * t35;
t19 = -qJDD(3) * pkin(3) - t82 * qJ(4) + t59 * t41 + qJDD(4) - t94;
t14 = (-t48 - t105) * pkin(7) + (t58 * t59 - qJDD(3)) * pkin(4) + t19;
t54 = -qJD(3) * pkin(4) - t59 * pkin(7);
t57 = t58 ^ 2;
t89 = -t82 * pkin(3) + qJDD(3) * qJ(4) + qJD(3) * t116 - t58 * t41 + t110;
t15 = -t57 * pkin(4) + t47 * pkin(7) + qJD(3) * t54 + t89;
t77 = sin(qJ(5));
t80 = cos(qJ(5));
t36 = t80 * t58 - t77 * t59;
t23 = t36 * qJD(5) + t77 * t47 + t80 * t48;
t37 = t77 * t58 + t80 * t59;
t26 = -t36 * mrSges(6,1) + t37 * mrSges(6,2);
t72 = -qJD(3) + qJD(5);
t30 = -t72 * mrSges(6,2) + t36 * mrSges(6,3);
t69 = -qJDD(3) + qJDD(5);
t10 = m(6) * (t80 * t14 - t77 * t15) - t23 * mrSges(6,3) + t69 * mrSges(6,1) - t37 * t26 + t72 * t30;
t22 = -t37 * qJD(5) + t80 * t47 - t77 * t48;
t31 = t72 * mrSges(6,1) - t37 * mrSges(6,3);
t11 = m(6) * (t77 * t14 + t80 * t15) + t22 * mrSges(6,3) - t69 * mrSges(6,2) + t36 * t26 - t72 * t31;
t52 = -qJD(3) * mrSges(5,1) + t59 * mrSges(5,2);
t91 = m(5) * t89 + qJDD(3) * mrSges(5,3) + qJD(3) * t52 - t77 * t10 + t80 * t11;
t6 = m(4) * t110 - qJDD(3) * mrSges(4,2) - qJD(3) * t51 + t109 * t58 + t111 * t47 + t91;
t50 = -qJD(3) * mrSges(4,2) - t58 * mrSges(4,3);
t53 = -t58 * mrSges(5,2) + qJD(3) * mrSges(5,3);
t88 = -m(5) * t19 - t80 * t10 - t77 * t11;
t7 = m(4) * t94 + t109 * t59 + t111 * t48 + (mrSges(4,1) + mrSges(5,1)) * qJDD(3) + (t50 + t53) * qJD(3) + t88;
t93 = -t76 * mrSges(3,1) + t75 * mrSges(3,2);
t92 = qJDD(1) * mrSges(3,3) + t83 * t93;
t4 = m(3) * t99 + t78 * t6 + t113 * t7 + (-m(3) * t60 - t92) * t75;
t5 = m(3) * t97 + t113 * t6 - t78 * t7 + t92 * t76;
t115 = t76 * t4 + t75 * t5;
t96 = m(6) * (-t57 * pkin(7) + (-pkin(3) - pkin(4)) * t47 + (-pkin(3) * qJD(3) + t116 + t54) * t59 - t122) + t23 * mrSges(6,2) - t22 * mrSges(6,1) + t37 * t31 - t36 * t30;
t87 = t48 * mrSges(5,3) + t59 * t52 - m(5) * (-0.2e1 * qJD(4) * t59 + (t47 + t104) * pkin(3) + t122) - t58 * t53 - t47 * mrSges(5,1) + t96;
t85 = m(4) * t86 + t47 * mrSges(4,1) + t48 * mrSges(4,2) + t58 * t50 + t59 * t51 - t87;
t84 = m(3) * (-qJDD(1) * pkin(1) - t83 * qJ(2) - t100) + t85;
t8 = -t84 + (-mrSges(2,2) + t121) * t83 + (mrSges(2,1) - t93) * qJDD(1) + m(2) * t108;
t1 = m(2) * t95 - t83 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t75 * t4 + t76 * t5;
t2 = [-m(1) * g(1) + t81 * t1 - t79 * t8, t1, t5, t6, -t47 * mrSges(5,2) - t58 * t42 + t91, t11; -m(1) * g(2) + t79 * t1 + t81 * t8, t8, t4, t7, -t87, t10; (-m(1) - m(2)) * g(3) + t115, -m(2) * g(3) + t115, t93 * qJDD(1) - t83 * t121 + t84, t85, -qJDD(3) * mrSges(5,1) + t48 * mrSges(5,2) - qJD(3) * t53 + t59 * t42 - t88, t96;];
f_new = t2;
