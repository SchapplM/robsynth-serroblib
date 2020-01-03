% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR7
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
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-31 18:20
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR7_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR7_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:18:58
% EndTime: 2019-12-31 18:19:00
% DurationCPUTime: 1.04s
% Computational Cost: add. (10879->138), mult. (23088->188), div. (0->0), fcn. (14369->10), ass. (0->74)
t65 = sin(pkin(9));
t67 = cos(pkin(9));
t70 = sin(qJ(3));
t73 = cos(qJ(3));
t45 = (t65 * t70 - t67 * t73) * qJD(1);
t92 = qJD(1) * qJD(3);
t88 = t73 * t92;
t54 = qJDD(1) * t70 + t88;
t55 = qJDD(1) * t73 - t70 * t92;
t94 = qJD(1) * t70;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t94;
t93 = qJD(1) * t73;
t58 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t93;
t76 = qJD(1) ^ 2;
t46 = (t65 * t73 + t67 * t70) * qJD(1);
t34 = pkin(4) * t45 - pkin(7) * t46;
t75 = qJD(3) ^ 2;
t71 = sin(qJ(1));
t74 = cos(qJ(1));
t89 = t71 * g(1) - g(2) * t74;
t51 = qJDD(1) * pkin(1) + t89;
t85 = -g(1) * t74 - g(2) * t71;
t53 = -pkin(1) * t76 + t85;
t66 = sin(pkin(8));
t68 = cos(pkin(8));
t95 = t66 * t51 + t68 * t53;
t31 = -pkin(2) * t76 + qJDD(1) * pkin(6) + t95;
t64 = -g(3) + qJDD(2);
t87 = -t70 * t31 + t73 * t64;
t20 = (-t54 + t88) * qJ(4) + (t70 * t73 * t76 + qJDD(3)) * pkin(3) + t87;
t56 = qJD(3) * pkin(3) - qJ(4) * t94;
t63 = t73 ^ 2;
t96 = t73 * t31 + t70 * t64;
t21 = -pkin(3) * t63 * t76 + qJ(4) * t55 - qJD(3) * t56 + t96;
t97 = 2 * qJD(4);
t90 = t65 * t20 + t67 * t21 - t45 * t97;
t16 = -pkin(4) * t75 + qJDD(3) * pkin(7) - t34 * t45 + t90;
t37 = -t54 * t65 + t55 * t67;
t38 = t54 * t67 + t55 * t65;
t86 = t68 * t51 - t66 * t53;
t81 = -qJDD(1) * pkin(2) - t86;
t78 = -t55 * pkin(3) + qJDD(4) + t56 * t94 + (-qJ(4) * t63 - pkin(6)) * t76 + t81;
t17 = (qJD(3) * t45 - t38) * pkin(7) + (qJD(3) * t46 - t37) * pkin(4) + t78;
t69 = sin(qJ(5));
t72 = cos(qJ(5));
t39 = qJD(3) * t72 - t46 * t69;
t25 = qJD(5) * t39 + qJDD(3) * t69 + t38 * t72;
t40 = qJD(3) * t69 + t46 * t72;
t26 = -mrSges(6,1) * t39 + mrSges(6,2) * t40;
t44 = qJD(5) + t45;
t27 = -mrSges(6,2) * t44 + mrSges(6,3) * t39;
t36 = qJDD(5) - t37;
t13 = m(6) * (-t16 * t69 + t17 * t72) - t25 * mrSges(6,3) + t36 * mrSges(6,1) - t40 * t26 + t44 * t27;
t24 = -qJD(5) * t40 + qJDD(3) * t72 - t38 * t69;
t28 = mrSges(6,1) * t44 - mrSges(6,3) * t40;
t14 = m(6) * (t16 * t72 + t17 * t69) + t24 * mrSges(6,3) - t36 * mrSges(6,2) + t39 * t26 - t44 * t28;
t41 = -qJD(3) * mrSges(5,2) - mrSges(5,3) * t45;
t42 = qJD(3) * mrSges(5,1) - mrSges(5,3) * t46;
t80 = -m(5) * t78 + t37 * mrSges(5,1) - t38 * mrSges(5,2) - t72 * t13 - t69 * t14 - t45 * t41 - t46 * t42;
t98 = (t57 * t70 - t58 * t73) * qJD(1) + m(4) * (-t76 * pkin(6) + t81) - t55 * mrSges(4,1) + t54 * mrSges(4,2) - t80;
t33 = mrSges(5,1) * t45 + mrSges(5,2) * t46;
t84 = -t67 * t20 + t65 * t21;
t79 = m(6) * (-qJDD(3) * pkin(4) - t75 * pkin(7) + (t97 + t34) * t46 + t84) - t24 * mrSges(6,1) + t25 * mrSges(6,2) - t39 * t27 + t40 * t28;
t10 = m(5) * (-0.2e1 * qJD(4) * t46 - t84) - t38 * mrSges(5,3) + qJDD(3) * mrSges(5,1) - t46 * t33 + qJD(3) * t41 - t79;
t52 = (-mrSges(4,1) * t73 + mrSges(4,2) * t70) * qJD(1);
t9 = m(5) * t90 - qJDD(3) * mrSges(5,2) + mrSges(5,3) * t37 - qJD(3) * t42 - t13 * t69 + t14 * t72 - t33 * t45;
t6 = m(4) * t87 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(3) * t58 + t67 * t10 - t52 * t94 + t65 * t9;
t7 = m(4) * t96 - qJDD(3) * mrSges(4,2) + mrSges(4,3) * t55 - qJD(3) * t57 - t10 * t65 + t52 * t93 + t67 * t9;
t91 = m(3) * t64 + t73 * t6 + t70 * t7;
t8 = m(3) * t86 + qJDD(1) * mrSges(3,1) - t76 * mrSges(3,2) - t98;
t3 = m(3) * t95 - mrSges(3,1) * t76 - qJDD(1) * mrSges(3,2) - t6 * t70 + t7 * t73;
t2 = m(2) * t85 - mrSges(2,1) * t76 - qJDD(1) * mrSges(2,2) + t3 * t68 - t66 * t8;
t1 = m(2) * t89 + qJDD(1) * mrSges(2,1) - mrSges(2,2) * t76 + t3 * t66 + t68 * t8;
t4 = [-m(1) * g(1) - t1 * t71 + t2 * t74, t2, t3, t7, t9, t14; -m(1) * g(2) + t1 * t74 + t2 * t71, t1, t8, t6, t10, t13; (-m(1) - m(2)) * g(3) + t91, -m(2) * g(3) + t91, t91, t98, -t80, t79;];
f_new = t4;
