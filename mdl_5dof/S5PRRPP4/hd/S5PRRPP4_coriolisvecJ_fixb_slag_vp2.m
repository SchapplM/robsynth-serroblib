% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,theta1]';
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:40:46
% EndTime: 2019-12-31 17:40:49
% DurationCPUTime: 1.00s
% Computational Cost: add. (524->193), mult. (1415->248), div. (0->0), fcn. (549->2), ass. (0->79)
t100 = Ifges(6,4) + Ifges(5,5);
t98 = mrSges(4,1) + mrSges(5,1);
t97 = -mrSges(5,2) - mrSges(4,3);
t85 = pkin(3) + pkin(4);
t55 = sin(qJ(3));
t71 = t55 * qJD(2);
t47 = pkin(6) * t71;
t56 = cos(qJ(3));
t31 = t56 * qJD(1) - t47;
t69 = qJ(5) * qJD(2);
t11 = t55 * t69 + t31;
t91 = -t11 + qJD(4);
t4 = -t85 * qJD(3) + t91;
t70 = t56 * qJD(2);
t32 = pkin(6) * t70 + t55 * qJD(1);
t12 = -t56 * t69 + t32;
t52 = qJD(3) * qJ(4);
t8 = t12 + t52;
t95 = m(6) * (t4 * t55 + t56 * t8);
t93 = -qJD(3) / 0.2e1;
t92 = qJD(3) / 0.2e1;
t90 = -t31 + qJD(4);
t89 = -t100 * t71 / 0.2e1;
t88 = 0.2e1 * pkin(6);
t87 = m(4) / 0.2e1;
t86 = m(5) / 0.2e1;
t84 = (pkin(2) * mrSges(4,1));
t83 = (pkin(2) * mrSges(4,2));
t82 = Ifges(4,4) * t55;
t27 = t32 * qJD(3);
t81 = t27 * t56;
t80 = mrSges(5,2) - mrSges(6,3);
t79 = pkin(6) - qJ(5);
t78 = t98 * qJD(3) + t97 * t71;
t74 = qJD(3) * mrSges(6,2);
t37 = -mrSges(6,3) * t70 + t74;
t39 = mrSges(5,2) * t70 + qJD(3) * mrSges(5,3);
t77 = t37 + t39;
t76 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t70 + t39;
t75 = qJ(4) * t56;
t73 = qJD(3) * t55;
t72 = qJD(4) * t55;
t68 = qJD(1) * qJD(3);
t41 = t79 * t56;
t66 = qJ(4) * t55 + pkin(2);
t65 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * Ifges(6,4);
t64 = -Ifges(6,5) / 0.2e1 + Ifges(5,4) / 0.2e1 + Ifges(4,5) / 0.2e1;
t63 = -Ifges(6,6) / 0.2e1 + Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1;
t61 = pkin(3) * t55 - t75;
t43 = t56 * t68;
t26 = -qJD(3) * t47 + t43;
t33 = -pkin(3) * t56 - t66;
t60 = -t85 * t55 + t75;
t24 = t85 * t56 + t66;
t15 = qJD(3) * t41 - qJD(5) * t55;
t13 = -qJD(5) * t56 - t79 * t73;
t14 = t61 * qJD(3) - t72;
t5 = t60 * qJD(3) + t72;
t16 = -qJD(3) * pkin(3) + t90;
t23 = t33 * qJD(2);
t46 = Ifges(4,4) * t70;
t6 = t24 * qJD(2) + qJD(5);
t59 = t16 * mrSges(5,2) + t6 * mrSges(6,2) + Ifges(6,5) * t93 + Ifges(4,1) * t71 / 0.2e1 + t46 / 0.2e1 - t23 * mrSges(5,3) - t31 * mrSges(4,3) - t4 * mrSges(6,3) + (-t100 * t56 + (Ifges(5,1) + Ifges(6,1)) * t55) * qJD(2) / 0.2e1 + (Ifges(5,4) + Ifges(4,5)) * t92;
t25 = t52 + t32;
t58 = t23 * mrSges(5,1) + t8 * mrSges(6,3) + Ifges(5,6) * t92 - (Ifges(4,2) * t56 + t82) * qJD(2) / 0.2e1 - t25 * mrSges(5,2) - t32 * mrSges(4,3) - t6 * mrSges(6,1) - t89 - (Ifges(5,3) + Ifges(6,2)) * t70 / 0.2e1 + (Ifges(6,6) + Ifges(4,6)) * t93;
t51 = qJD(3) * qJD(4);
t42 = t70 * t74;
t40 = t79 * t55;
t34 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t71;
t30 = t61 * qJD(2);
t29 = (t56 * mrSges(6,1) + t55 * mrSges(6,2)) * qJD(2);
t28 = (-mrSges(5,1) * t56 - mrSges(5,3) * t55) * qJD(2);
t10 = t26 + t51;
t9 = t60 * qJD(2);
t7 = t14 * qJD(2);
t3 = t15 * qJD(2) + t55 * t68;
t2 = t5 * qJD(2);
t1 = t13 * qJD(2) + t43 + t51;
t17 = [m(4) * (t26 * t55 - t81) + m(5) * (t10 * t55 - t81) + m(6) * (t1 * t55 - t3 * t56) + ((t37 + t76) * t56 + (t34 - t78) * t55 + m(4) * (-t31 * t55 + t32 * t56) + m(5) * (t16 * t55 + t25 * t56) + t95 + (t55 ^ 2 + t56 ^ 2) * qJD(2) * (-mrSges(4,3) - t80)) * qJD(3); t13 * t37 + t14 * t28 + t15 * t34 + t24 * t42 + t5 * t29 + m(5) * (t14 * t23 + t33 * t7) + m(6) * (t1 * t41 + t13 * t8 + t15 * t4 + t2 * t24 + t3 * t40 + t5 * t6) + (-t7 * mrSges(5,1) + t2 * mrSges(6,1) + t10 * mrSges(5,2) + t26 * mrSges(4,3) - t1 * mrSges(6,3) + (t10 * t86 + t26 * t87) * t88 + (t64 * qJD(3) + (-m(4) * t31 + m(5) * t16 - t78) * pkin(6) + (-t33 * mrSges(5,3) - t40 * mrSges(6,3) - t65 * t56 - (2 * t83)) * qJD(2) + t59) * qJD(3)) * t56 + (t2 * mrSges(6,2) - t7 * mrSges(5,3) - t3 * mrSges(6,3) + ((t87 + t86) * t88 - t97) * t27 + (t63 * qJD(3) + (-m(4) * t32 - m(5) * t25 - t76) * pkin(6) + (t41 * mrSges(6,3) + t33 * mrSges(5,1) - (2 * t84) - t24 * mrSges(6,1) + t65 * t55 + (-0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2) - 0.3e1 / 0.2e1 * Ifges(6,2) + 0.3e1 / 0.2e1 * Ifges(6,1) + 0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1)) * t56) * qJD(2) + t58) * qJD(3)) * t55; -t3 * mrSges(6,1) - t26 * mrSges(4,2) + t1 * mrSges(6,2) + t10 * mrSges(5,3) - t11 * t37 - t12 * t34 - t30 * t28 - t9 * t29 + t78 * t32 - t76 * t31 - t98 * t27 + t77 * qJD(4) + (((t84 + t82 / 0.2e1) * qJD(2) + (-t80 * qJ(4) + t63) * qJD(3) - t58 + t89) * t55 + (-t46 / 0.2e1 + (t83 + (Ifges(5,5) / 0.2e1 + Ifges(6,4) / 0.2e1) * t56) * qJD(2) + (Ifges(5,3) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1 - Ifges(5,1) / 0.2e1 - Ifges(6,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t71 + (-pkin(3) * mrSges(5,2) + mrSges(6,3) * t85 + t64) * qJD(3) - t59) * t56) * qJD(2) + (t1 * qJ(4) - t4 * t12 - t3 * t85 - t6 * t9 + t91 * t8) * m(6) + (-t27 * pkin(3) + t10 * qJ(4) - t16 * t32 - t23 * t30 + t90 * t25) * m(5); -t77 * qJD(3) + ((t28 - t29) * t55 + t80 * t56 * qJD(3)) * qJD(2) + (-t8 * qJD(3) - t6 * t71 + t3) * m(6) + (-t25 * qJD(3) + t23 * t71 + t27) * m(5); t42 + m(6) * t2 + (-mrSges(6,1) * t73 + t55 * t34 + t56 * t37 + t95) * qJD(2);];
tauc = t17(:);
