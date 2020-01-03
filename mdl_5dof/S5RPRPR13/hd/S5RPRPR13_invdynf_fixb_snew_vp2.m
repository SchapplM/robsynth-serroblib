% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR13
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR13_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:59
% EndTime: 2019-12-31 18:32:02
% DurationCPUTime: 0.94s
% Computational Cost: add. (7258->151), mult. (17417->186), div. (0->0), fcn. (11702->8), ass. (0->79)
t77 = qJD(1) ^ 2;
t108 = cos(qJ(3));
t69 = sin(pkin(8));
t70 = cos(pkin(8));
t72 = sin(qJ(3));
t113 = -t70 * t108 + t69 * t72;
t58 = t113 * qJD(1);
t49 = t58 * mrSges(5,1) - qJD(3) * mrSges(5,3);
t102 = -qJD(3) * mrSges(4,2) - t58 * mrSges(4,3) - t49;
t106 = mrSges(4,1) - mrSges(5,2);
t86 = t108 * t69 + t70 * t72;
t59 = t86 * qJD(1);
t98 = t59 * qJD(3);
t43 = t113 * qJDD(1) + t98;
t99 = t58 * qJD(3);
t44 = t86 * qJDD(1) - t99;
t48 = qJD(3) * mrSges(4,1) - t59 * mrSges(4,3);
t68 = t70 ^ 2;
t101 = t69 ^ 2 + t68;
t73 = sin(qJ(1));
t75 = cos(qJ(1));
t95 = t73 * g(1) - t75 * g(2);
t91 = qJDD(2) - t95;
t79 = (-pkin(2) * t70 - pkin(1)) * qJDD(1) + (-t101 * pkin(6) - qJ(2)) * t77 + t91;
t51 = t59 * pkin(4) - qJD(3) * pkin(7);
t57 = t58 ^ 2;
t112 = -2 * qJD(4);
t78 = pkin(3) * t98 + t59 * t112 + (-t44 + t99) * qJ(4) + t79;
t12 = -t57 * pkin(4) - t59 * t51 + (pkin(3) + pkin(7)) * t43 + t78;
t36 = t58 * pkin(3) - t59 * qJ(4);
t76 = qJD(3) ^ 2;
t100 = pkin(6) * qJDD(1);
t109 = pkin(2) * t77;
t92 = -t75 * g(1) - t73 * g(2);
t60 = -t77 * pkin(1) + qJDD(1) * qJ(2) + t92;
t97 = qJD(1) * qJD(2);
t94 = -t70 * g(3) - 0.2e1 * t69 * t97;
t31 = (t70 * t109 - t100 - t60) * t69 + t94;
t93 = -t69 * g(3) + (t60 + 0.2e1 * t97) * t70;
t32 = t70 * t100 - t68 * t109 + t93;
t90 = t108 * t31 - t72 * t32;
t19 = -qJDD(3) * pkin(3) - t76 * qJ(4) + t59 * t36 + qJDD(4) - t90;
t13 = (t58 * t59 - qJDD(3)) * pkin(7) + (t44 + t99) * pkin(4) + t19;
t71 = sin(qJ(5));
t74 = cos(qJ(5));
t45 = -t71 * qJD(3) + t74 * t58;
t23 = t45 * qJD(5) + t74 * qJDD(3) + t71 * t43;
t46 = t74 * qJD(3) + t71 * t58;
t24 = -t45 * mrSges(6,1) + t46 * mrSges(6,2);
t55 = qJD(5) + t59;
t27 = -t55 * mrSges(6,2) + t45 * mrSges(6,3);
t41 = qJDD(5) + t44;
t10 = m(6) * (-t71 * t12 + t74 * t13) - t23 * mrSges(6,3) + t41 * mrSges(6,1) - t46 * t24 + t55 * t27;
t22 = -t46 * qJD(5) - t71 * qJDD(3) + t74 * t43;
t28 = t55 * mrSges(6,1) - t46 * mrSges(6,3);
t11 = m(6) * (t74 * t12 + t71 * t13) + t22 * mrSges(6,3) - t41 * mrSges(6,2) + t45 * t24 - t55 * t28;
t50 = t59 * mrSges(5,1) + qJD(3) * mrSges(5,2);
t87 = t71 * t10 - t74 * t11 - m(5) * (t43 * pkin(3) + t78) + t59 * t50 + t44 * mrSges(5,3);
t80 = m(4) * t79 + t44 * mrSges(4,2) + t102 * t58 + t106 * t43 + t59 * t48 - t87;
t115 = -m(3) * (-qJDD(1) * pkin(1) - t77 * qJ(2) + t91) - t80;
t114 = t101 * mrSges(3,3);
t38 = -t58 * mrSges(5,2) - t59 * mrSges(5,3);
t103 = -t58 * mrSges(4,1) - t59 * mrSges(4,2) - t38;
t105 = -mrSges(4,3) - mrSges(5,1);
t85 = -m(5) * t19 - t74 * t10 - t71 * t11;
t7 = m(4) * t90 + t102 * qJD(3) + t106 * qJDD(3) + t103 * t59 + t105 * t44 + t85;
t104 = t108 * t32 + t72 * t31;
t82 = -t76 * pkin(3) + qJDD(3) * qJ(4) - t58 * t36 + t104;
t84 = -t22 * mrSges(6,1) - t45 * t27 + m(6) * (-t43 * pkin(4) - t57 * pkin(7) + ((2 * qJD(4)) + t51) * qJD(3) + t82) + t23 * mrSges(6,2) + t46 * t28;
t81 = -m(5) * (qJD(3) * t112 - t82) + t84;
t8 = m(4) * t104 + t103 * t58 + t105 * t43 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) + (-t48 + t50) * qJD(3) + t81;
t89 = -t70 * mrSges(3,1) + t69 * mrSges(3,2);
t88 = qJDD(1) * mrSges(3,3) + t77 * t89;
t4 = m(3) * t94 + t72 * t8 + t108 * t7 + (-m(3) * t60 - t88) * t69;
t5 = m(3) * t93 + t108 * t8 - t72 * t7 + t88 * t70;
t111 = t70 * t4 + t69 * t5;
t6 = m(2) * t95 + (-mrSges(2,2) + t114) * t77 + (mrSges(2,1) - t89) * qJDD(1) + t115;
t1 = m(2) * t92 - t77 * mrSges(2,1) - qJDD(1) * mrSges(2,2) - t69 * t4 + t70 * t5;
t2 = [-m(1) * g(1) + t75 * t1 - t73 * t6, t1, t5, t8, -t43 * mrSges(5,2) - t58 * t49 - t87, t11; -m(1) * g(2) + t73 * t1 + t75 * t6, t6, t4, t7, t43 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t50 + t58 * t38 - t81, t10; (-m(1) - m(2)) * g(3) + t111, -m(2) * g(3) + t111, t89 * qJDD(1) - t77 * t114 - t115, t80, t44 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t49 + t59 * t38 - t85, t84;];
f_new = t2;
