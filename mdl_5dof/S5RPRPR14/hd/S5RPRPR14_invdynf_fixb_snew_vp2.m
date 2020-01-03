% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR14
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
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta4]';
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
% Datum: 2019-12-31 18:35
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR14_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR14_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:34:33
% EndTime: 2019-12-31 18:34:35
% DurationCPUTime: 0.80s
% Computational Cost: add. (7156->138), mult. (15342->180), div. (0->0), fcn. (9352->8), ass. (0->72)
t68 = sin(qJ(1));
t71 = cos(qJ(1));
t85 = -t71 * g(1) - t68 * g(2);
t81 = qJDD(1) * qJ(2) + (2 * qJD(2) * qJD(1)) + t85;
t64 = sin(pkin(8));
t65 = cos(pkin(8));
t67 = sin(qJ(3));
t70 = cos(qJ(3));
t46 = (t64 * t70 + t65 * t67) * qJD(1);
t100 = 2 * qJD(4);
t99 = -m(2) - m(3);
t98 = -pkin(1) - pkin(6);
t97 = (mrSges(2,1) - mrSges(3,2));
t96 = -mrSges(2,2) + mrSges(3,3);
t73 = qJD(1) ^ 2;
t88 = t68 * g(1) - t71 * g(2);
t79 = -t73 * qJ(2) + qJDD(2) - t88;
t41 = t98 * qJDD(1) + t79;
t95 = t67 * g(3) + t70 * t41;
t94 = qJD(1) * t67;
t93 = qJD(1) * t70;
t92 = qJD(1) * qJD(3);
t86 = t67 * t92;
t54 = t70 * qJDD(1) - t86;
t19 = (-t54 - t86) * qJ(4) + (-t67 * t70 * t73 + qJDD(3)) * pkin(3) + t95;
t53 = -t67 * qJDD(1) - t70 * t92;
t56 = qJD(3) * pkin(3) - qJ(4) * t93;
t63 = t67 ^ 2;
t87 = -t70 * g(3) + t67 * t41;
t20 = -t63 * t73 * pkin(3) + t53 * qJ(4) - qJD(3) * t56 + t87;
t90 = -t46 * t100 + t64 * t19 + t65 * t20;
t52 = (mrSges(4,1) * t67 + mrSges(4,2) * t70) * qJD(1);
t55 = -qJD(3) * mrSges(4,2) - mrSges(4,3) * t94;
t47 = -t64 * t94 + t65 * t93;
t28 = t46 * pkin(4) - t47 * pkin(7);
t72 = qJD(3) ^ 2;
t13 = -t72 * pkin(4) + qJDD(3) * pkin(7) - t46 * t28 + t90;
t31 = t65 * t53 - t64 * t54;
t32 = t64 * t53 + t65 * t54;
t75 = -t53 * pkin(3) + qJDD(4) + t56 * t93 + (-qJ(4) * t63 + t98) * t73 + t81;
t14 = (qJD(3) * t46 - t32) * pkin(7) + (qJD(3) * t47 - t31) * pkin(4) + t75;
t66 = sin(qJ(5));
t69 = cos(qJ(5));
t33 = t69 * qJD(3) - t66 * t47;
t16 = t33 * qJD(5) + t66 * qJDD(3) + t69 * t32;
t34 = t66 * qJD(3) + t69 * t47;
t23 = -t33 * mrSges(6,1) + t34 * mrSges(6,2);
t44 = qJD(5) + t46;
t24 = -t44 * mrSges(6,2) + t33 * mrSges(6,3);
t30 = qJDD(5) - t31;
t10 = m(6) * (-t66 * t13 + t69 * t14) - t16 * mrSges(6,3) + t30 * mrSges(6,1) - t34 * t23 + t44 * t24;
t15 = -t34 * qJD(5) + t69 * qJDD(3) - t66 * t32;
t25 = t44 * mrSges(6,1) - t34 * mrSges(6,3);
t11 = m(6) * (t69 * t13 + t66 * t14) + t15 * mrSges(6,3) - t30 * mrSges(6,2) + t33 * t23 - t44 * t25;
t27 = t46 * mrSges(5,1) + t47 * mrSges(5,2);
t40 = qJD(3) * mrSges(5,1) - t47 * mrSges(5,3);
t6 = m(5) * t90 - qJDD(3) * mrSges(5,2) + t31 * mrSges(5,3) - qJD(3) * t40 - t66 * t10 + t69 * t11 - t46 * t27;
t39 = -qJD(3) * mrSges(5,2) - t46 * mrSges(5,3);
t84 = -t65 * t19 + t64 * t20;
t77 = m(6) * (-qJDD(3) * pkin(4) - t72 * pkin(7) + (t100 + t28) * t47 + t84) - t15 * mrSges(6,1) + t16 * mrSges(6,2) - t33 * t24 + t34 * t25;
t7 = m(5) * (-0.2e1 * qJD(4) * t47 - t84) - t32 * mrSges(5,3) + qJDD(3) * mrSges(5,1) - t47 * t27 + qJD(3) * t39 - t77;
t3 = m(4) * t95 + qJDD(3) * mrSges(4,1) - t54 * mrSges(4,3) + qJD(3) * t55 - t52 * t93 + t64 * t6 + t65 * t7;
t57 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t93;
t4 = m(4) * t87 - qJDD(3) * mrSges(4,2) + t53 * mrSges(4,3) - qJD(3) * t57 - t52 * t94 + t65 * t6 - t64 * t7;
t89 = -t67 * t3 + t70 * t4;
t80 = -m(3) * (-qJDD(1) * pkin(1) + t79) - t70 * t3 - t67 * t4;
t78 = m(5) * t75 - t31 * mrSges(5,1) + t32 * mrSges(5,2) + t69 * t10 + t66 * t11 + t46 * t39 + t47 * t40;
t76 = -t53 * mrSges(4,1) + m(4) * (t98 * t73 + t81) + t55 * t94 + t57 * t93 + t54 * mrSges(4,2) + t78;
t74 = -m(3) * (t73 * pkin(1) - t81) + t76;
t5 = m(2) * t85 + t96 * qJDD(1) - (t97 * t73) + t74;
t1 = m(2) * t88 + t97 * qJDD(1) + t96 * t73 + t80;
t2 = [-m(1) * g(1) - t68 * t1 + t71 * t5, t5, -m(3) * g(3) + t89, t4, t6, t11; -m(1) * g(2) + t71 * t1 + t68 * t5, t1, -(t73 * mrSges(3,2)) - qJDD(1) * mrSges(3,3) - t74, t3, t7, t10; (-m(1) + t99) * g(3) + t89, t99 * g(3) + t89, qJDD(1) * mrSges(3,2) - t73 * mrSges(3,3) - t80, t76, t78, t77;];
f_new = t2;
