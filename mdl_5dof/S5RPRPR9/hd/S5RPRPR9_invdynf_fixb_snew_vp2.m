% Calculate vector of cutting forces with Newton-Euler
% S5RPRPR9
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
% Datum: 2019-12-31 18:25
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function f_new = S5RPRPR9_invdynf_fixb_snew_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR9_invdynf_fixb_snew_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_NewtonEuler_linkframe_f_i_i_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:23:45
% EndTime: 2019-12-31 18:23:46
% DurationCPUTime: 0.64s
% Computational Cost: add. (4438->139), mult. (8634->172), div. (0->0), fcn. (4346->8), ass. (0->73)
t64 = sin(qJ(3));
t67 = cos(qJ(3));
t86 = qJD(1) * qJD(3);
t82 = t67 * t86;
t45 = t64 * qJDD(1) + t82;
t83 = t64 * t86;
t46 = t67 * qJDD(1) - t83;
t88 = qJD(1) * t67;
t48 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t88;
t65 = sin(qJ(1));
t68 = cos(qJ(1));
t84 = t65 * g(1) - t68 * g(2);
t40 = qJDD(1) * pkin(1) + t84;
t70 = qJD(1) ^ 2;
t80 = -t68 * g(1) - t65 * g(2);
t44 = -t70 * pkin(1) + t80;
t61 = sin(pkin(8));
t62 = cos(pkin(8));
t81 = t62 * t40 - t61 * t44;
t77 = -qJDD(1) * pkin(2) - t81;
t87 = t64 * qJD(1);
t51 = pkin(4) * t87 - qJD(3) * pkin(7);
t59 = t67 ^ 2;
t100 = -2 * qJD(4);
t71 = pkin(3) * t83 + t87 * t100 + (-t45 - t82) * qJ(4) + t77;
t99 = -pkin(3) - pkin(7);
t12 = -t51 * t87 + (-pkin(4) * t59 - pkin(6)) * t70 + t99 * t46 + t71;
t60 = -g(3) + qJDD(2);
t91 = t61 * t40 + t62 * t44;
t23 = -t70 * pkin(2) + qJDD(1) * pkin(6) + t91;
t20 = t64 * t23;
t41 = (-pkin(3) * t67 - qJ(4) * t64) * qJD(1);
t69 = qJD(3) ^ 2;
t79 = -t69 * qJ(4) + t41 * t87 + qJDD(4) + t20;
t98 = pkin(7) * t70;
t15 = t45 * pkin(4) + t99 * qJDD(3) + (-pkin(4) * t86 - t64 * t98 - t60) * t67 + t79;
t63 = sin(qJ(5));
t66 = cos(qJ(5));
t38 = -t63 * qJD(3) - t66 * t88;
t27 = t38 * qJD(5) + t66 * qJDD(3) - t63 * t46;
t39 = t66 * qJD(3) - t63 * t88;
t28 = -t38 * mrSges(6,1) + t39 * mrSges(6,2);
t53 = qJD(5) + t87;
t29 = -t53 * mrSges(6,2) + t38 * mrSges(6,3);
t37 = qJDD(5) + t45;
t10 = m(6) * (-t63 * t12 + t66 * t15) - t27 * mrSges(6,3) + t37 * mrSges(6,1) - t39 * t28 + t53 * t29;
t26 = -t39 * qJD(5) - t63 * qJDD(3) - t66 * t46;
t30 = t53 * mrSges(6,1) - t39 * mrSges(6,3);
t11 = m(6) * (t66 * t12 + t63 * t15) + t26 * mrSges(6,3) - t37 * mrSges(6,2) + t38 * t28 - t53 * t30;
t49 = -mrSges(5,1) * t88 - qJD(3) * mrSges(5,3);
t97 = t70 * pkin(6);
t78 = t63 * t10 - t66 * t11 - m(5) * (-t46 * pkin(3) + t71 - t97) - t49 * t88 + t45 * mrSges(5,3);
t50 = mrSges(5,1) * t87 + qJD(3) * mrSges(5,2);
t89 = qJD(3) * mrSges(4,1) - mrSges(4,3) * t87 - t50;
t94 = mrSges(4,1) - mrSges(5,2);
t103 = (-t67 * t48 + t89 * t64) * qJD(1) - t94 * t46 + m(4) * (t77 - t97) + t45 * mrSges(4,2) - t78;
t95 = t67 * t60;
t93 = mrSges(4,3) + mrSges(5,1);
t92 = t67 * t23 + t64 * t60;
t42 = (mrSges(5,2) * t67 - mrSges(5,3) * t64) * qJD(1);
t90 = t42 + (-mrSges(4,1) * t67 + mrSges(4,2) * t64) * qJD(1);
t76 = -m(5) * (-qJDD(3) * pkin(3) + t79 - t95) - t66 * t10 - t63 * t11;
t6 = m(4) * (-t20 + t95) - t93 * t45 + t94 * qJDD(3) + (t48 - t49) * qJD(3) - t90 * t87 + t76;
t72 = -t69 * pkin(3) + qJDD(3) * qJ(4) + t41 * t88 + t92;
t75 = -t26 * mrSges(6,1) - t38 * t29 + m(6) * (-t59 * t98 + t46 * pkin(4) + ((2 * qJD(4)) + t51) * qJD(3) + t72) + t27 * mrSges(6,2) + t39 * t30;
t74 = -m(5) * (qJD(3) * t100 - t72) + t75;
t8 = m(4) * t92 + t93 * t46 + (-mrSges(4,2) + mrSges(5,3)) * qJDD(3) - t89 * qJD(3) + t90 * t88 + t74;
t85 = m(3) * t60 + t67 * t6 + t64 * t8;
t4 = m(3) * t81 + qJDD(1) * mrSges(3,1) - t70 * mrSges(3,2) - t103;
t3 = m(3) * t91 - t70 * mrSges(3,1) - qJDD(1) * mrSges(3,2) - t64 * t6 + t67 * t8;
t2 = m(2) * t80 - t70 * mrSges(2,1) - qJDD(1) * mrSges(2,2) + t62 * t3 - t61 * t4;
t1 = m(2) * t84 + qJDD(1) * mrSges(2,1) - t70 * mrSges(2,2) + t61 * t3 + t62 * t4;
t5 = [-m(1) * g(1) - t65 * t1 + t68 * t2, t2, t3, t8, t46 * mrSges(5,2) - t50 * t87 - t78, t11; -m(1) * g(2) + t68 * t1 + t65 * t2, t1, t4, t6, -t46 * mrSges(5,1) - qJDD(3) * mrSges(5,3) - qJD(3) * t50 - t42 * t88 - t74, t10; (-m(1) - m(2)) * g(3) + t85, -m(2) * g(3) + t85, t85, t103, t45 * mrSges(5,1) + qJDD(3) * mrSges(5,2) + qJD(3) * t49 + t42 * t87 - t76, t75;];
f_new = t5;
