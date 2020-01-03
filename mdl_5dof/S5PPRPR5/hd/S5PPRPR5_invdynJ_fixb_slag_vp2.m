% Calculate vector of inverse dynamics joint torques for
% S5PPRPR5
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
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% tau [5x1]
%   joint torques of inverse dynamics (contains inertial, gravitational coriolis and centrifugal forces)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tau = S5PPRPR5_invdynJ_fixb_slag_vp2(qJ, qJD, qJDD, g, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(3,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(qJDD) && all(size(qJDD) == [5 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: qJDD has to be [5x1] (double)');
assert(isreal(g) && all(size(g) == [3 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: g has to be [3x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_invdynJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From invdyn_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:25
% EndTime: 2019-12-31 17:33:26
% DurationCPUTime: 1.08s
% Computational Cost: add. (536->157), mult. (1012->197), div. (0->0), fcn. (528->6), ass. (0->79)
t107 = -m(6) - m(5);
t32 = sin(qJ(5));
t106 = -t32 / 0.2e1;
t34 = cos(qJ(5));
t91 = t34 / 0.2e1;
t33 = sin(qJ(3));
t35 = cos(qJ(3));
t74 = sin(pkin(7));
t75 = cos(pkin(7));
t16 = -t33 * t74 - t35 * t75;
t17 = t33 * t75 - t35 * t74;
t99 = -g(1) * t17 + g(2) * t16;
t36 = -pkin(3) - pkin(6);
t72 = qJD(2) * t35;
t56 = qJD(4) - t72;
t39 = qJD(3) * t36 + t56;
t4 = t32 * qJD(1) + t34 * t39;
t5 = t34 * qJD(1) - t32 * t39;
t63 = qJD(2) * qJD(3);
t31 = t33 * t63;
t21 = qJDD(2) * t35 - t31;
t46 = qJDD(4) - t21;
t6 = qJDD(3) * t36 + t46;
t1 = qJD(5) * t4 - t34 * qJDD(1) + t32 * t6;
t2 = qJD(5) * t5 + t32 * qJDD(1) + t34 * t6;
t55 = t1 * t32 + t2 * t34;
t69 = qJD(5) * t34;
t70 = qJD(5) * t32;
t105 = t4 * t70 + t5 * t69 - t55;
t77 = mrSges(4,1) - mrSges(5,2);
t104 = -m(6) * pkin(6) + pkin(3) * t107 - t77;
t71 = qJD(3) * t32;
t24 = -qJD(5) * mrSges(6,2) - mrSges(6,3) * t71;
t66 = t34 * qJD(3);
t25 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t66;
t101 = t32 * t24 + t34 * t25;
t48 = t34 * t24 - t32 * t25;
t62 = qJD(3) * qJD(5);
t19 = qJDD(3) * t34 - t32 * t62;
t8 = qJDD(5) * mrSges(6,1) - mrSges(6,3) * t19;
t20 = -qJDD(3) * t32 - t34 * t62;
t9 = -qJDD(5) * mrSges(6,2) + mrSges(6,3) * t20;
t100 = t32 * t9 + t34 * t8;
t98 = m(3) + m(4) + m(5);
t83 = Ifges(6,4) * t34;
t84 = Ifges(6,4) * t32;
t97 = (-Ifges(6,2) * t34 - t84) * t106 + (-Ifges(6,1) * t32 - t83) * t91;
t26 = t32 * mrSges(6,1) + t34 * mrSges(6,2);
t18 = t26 * qJD(3);
t95 = -qJD(3) * t18 + t48 * qJD(5) + t100;
t73 = qJD(2) * t33;
t27 = qJD(3) * qJ(4) + t73;
t52 = mrSges(6,1) * t34 - mrSges(6,2) * t32;
t93 = t27 * t52 + qJD(5) * (-Ifges(6,5) * t32 - Ifges(6,6) * t34) / 0.2e1;
t92 = -t107 * qJ(4) - mrSges(4,2) + t26;
t82 = t27 * t35;
t76 = -mrSges(4,2) + mrSges(5,3);
t67 = t27 * qJD(3);
t65 = qJDD(2) * t33;
t64 = qJDD(3) * mrSges(5,2);
t61 = qJDD(3) * qJ(4);
t58 = m(2) + t98;
t57 = t35 * t63;
t54 = -t32 * t5 + t34 * t4;
t53 = -t4 * t32 - t5 * t34;
t51 = Ifges(6,1) * t34 - t84;
t50 = -Ifges(6,2) * t32 + t83;
t7 = t61 + t65 + (qJD(4) + t72) * qJD(3);
t47 = qJ(4) * t7 + t27 * qJD(4);
t41 = -t67 + t99;
t38 = qJD(5) * t53 + t55;
t37 = qJD(3) ^ 2;
t23 = -qJD(3) * pkin(3) + t56;
t22 = t57 + t65;
t13 = Ifges(6,5) * qJD(5) + qJD(3) * t51;
t12 = Ifges(6,6) * qJD(5) + qJD(3) * t50;
t10 = -qJDD(3) * pkin(3) + t46;
t3 = -mrSges(6,1) * t20 + mrSges(6,2) * t19;
t11 = [t24 * t70 + t25 * t69 - t34 * t9 + t32 * t8 + m(6) * (qJD(5) * t54 - t1 * t34 + t2 * t32) + t58 * qJDD(1) + (-m(6) - t58) * g(3); m(3) * qJDD(2) + (t3 - t77 * t37 + t76 * qJDD(3) + t101 * qJD(3) + m(4) * t22 + m(5) * (qJD(3) * t23 + t7) + m(6) * (t4 * t66 - t5 * t71 + t7)) * t33 + (t76 * t37 + t77 * qJDD(3) + m(4) * t21 + m(5) * (-t10 + t67) + m(6) * (t67 + t105) - t95) * t35 + (-t74 * g(1) + t75 * g(2)) * (m(6) + t98); (Ifges(5,1) + Ifges(4,3)) * qJDD(3) + (t105 - t99) * mrSges(6,3) + (0.2e1 * Ifges(6,5) * t91 - Ifges(6,6) * t32) * qJDD(5) + t56 * t18 + t100 * t36 - t101 * t73 + t97 * t62 + (t104 * t16 + t92 * t17) * g(2) + (-t104 * t17 + t92 * t16) * g(1) + (t48 * t36 + t93) * qJD(5) + (Ifges(6,1) * t19 + Ifges(6,4) * t20) * t91 + (g(1) * t16 + g(2) * t17 + qJD(3) * qJD(4) - t57 + t61 + t7) * mrSges(5,3) + t7 * t26 + qJ(4) * t3 + (t10 - t31) * mrSges(5,2) + (-t22 + t57) * mrSges(4,2) + (-(t33 * t54 + t82) * qJD(2) + t38 * t36 + t47) * m(6) + (-(t23 * t33 + t82) * qJD(2) - pkin(3) * t10 + t47) * m(5) + (Ifges(6,4) * t19 + Ifges(6,2) * t20) * t106 + (t21 + t31) * mrSges(4,1) + t20 * t50 / 0.2e1 + t19 * t51 / 0.2e1 - pkin(3) * t64 - t12 * t69 / 0.2e1 - t13 * t70 / 0.2e1; t64 - t37 * mrSges(5,3) + (t38 + t41) * m(6) + (t10 + t41) * m(5) + t95; t2 * mrSges(6,1) - t1 * mrSges(6,2) + Ifges(6,5) * t19 + Ifges(6,6) * t20 + Ifges(6,3) * qJDD(5) - g(3) * t26 - t4 * t24 - t5 * t25 + (t32 * t13 / 0.2e1 + t12 * t91 - t97 * qJD(3) + t53 * mrSges(6,3) - t93) * qJD(3) + t99 * t52;];
tau = t11;
