% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,theta2]';
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
% Datum: 2019-12-31 18:11
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:10:41
% EndTime: 2019-12-31 18:10:43
% DurationCPUTime: 1.03s
% Computational Cost: add. (639->197), mult. (1620->251), div. (0->0), fcn. (664->4), ass. (0->80)
t104 = Ifges(6,4) + Ifges(5,5);
t102 = -mrSges(4,1) - mrSges(5,1);
t101 = mrSges(5,2) + mrSges(4,3);
t60 = sin(qJ(3));
t76 = t60 * qJD(1);
t99 = t76 / 0.2e1;
t88 = pkin(3) + pkin(4);
t47 = sin(pkin(7)) * pkin(1) + pkin(6);
t38 = t47 * qJD(1);
t61 = cos(qJ(3));
t53 = t61 * qJD(2);
t20 = -t60 * t38 + t53;
t74 = qJ(5) * qJD(1);
t8 = t60 * t74 + t20;
t94 = -t8 + qJD(4);
t3 = -t88 * qJD(3) + t94;
t55 = qJD(3) * qJ(4);
t52 = t60 * qJD(2);
t21 = t61 * t38 + t52;
t9 = -t61 * t74 + t21;
t6 = t55 + t9;
t98 = m(6) * (t3 * t60 + t6 * t61);
t96 = -qJD(3) / 0.2e1;
t95 = qJD(3) / 0.2e1;
t93 = -t20 + qJD(4);
t92 = -t104 * t76 / 0.2e1;
t91 = 0.2e1 * t47;
t90 = m(4) / 0.2e1;
t89 = m(5) / 0.2e1;
t18 = t21 * qJD(3);
t87 = t18 * t61;
t86 = mrSges(5,2) - mrSges(6,3);
t85 = t102 * qJD(3) + t101 * t76;
t75 = t61 * qJD(1);
t80 = qJD(3) * mrSges(6,2);
t43 = -mrSges(6,3) * t75 + t80;
t45 = mrSges(5,2) * t75 + qJD(3) * mrSges(5,3);
t84 = t43 + t45;
t83 = -qJD(3) * mrSges(4,2) + mrSges(4,3) * t75 + t45;
t82 = qJ(4) * t61;
t81 = qJ(5) - t47;
t79 = qJD(3) * t60;
t78 = qJD(4) * t60;
t77 = qJD(5) * t61;
t72 = -cos(pkin(7)) * pkin(1) - pkin(2);
t32 = t81 * t61;
t71 = Ifges(4,5) / 0.2e1 - Ifges(6,5) / 0.2e1 + Ifges(5,4) / 0.2e1;
t70 = 0.3e1 / 0.2e1 * Ifges(5,5) - 0.3e1 / 0.2e1 * Ifges(4,4) + 0.3e1 / 0.2e1 * Ifges(6,4);
t69 = Ifges(5,6) / 0.2e1 - Ifges(4,6) / 0.2e1 - Ifges(6,6) / 0.2e1;
t17 = qJD(3) * t53 - t38 * t79;
t67 = pkin(3) * t60 - t82;
t7 = qJD(3) * qJD(4) + t17;
t66 = qJ(4) * t60 - t72;
t65 = -t88 * t60 + t82;
t30 = -pkin(3) * t61 - t66;
t19 = t88 * t61 + t66;
t23 = qJD(3) * t67 - t78;
t11 = qJD(3) * t65 + t78;
t10 = -qJD(3) * pkin(3) + t93;
t15 = t30 * qJD(1);
t4 = qJD(1) * t19 + qJD(5);
t42 = t72 * qJD(1);
t51 = Ifges(4,4) * t75;
t64 = t10 * mrSges(5,2) + t4 * mrSges(6,2) + t42 * mrSges(4,2) + Ifges(6,5) * t96 + Ifges(4,1) * t99 + t51 / 0.2e1 - t15 * mrSges(5,3) - t20 * mrSges(4,3) - t3 * mrSges(6,3) + (-t104 * t61 + (Ifges(5,1) + Ifges(6,1)) * t60) * qJD(1) / 0.2e1 + (Ifges(5,4) + Ifges(4,5)) * t95;
t14 = t55 + t21;
t63 = t15 * mrSges(5,1) + t42 * mrSges(4,1) + t6 * mrSges(6,3) + Ifges(5,6) * t95 - (Ifges(4,4) * t60 + Ifges(4,2) * t61) * qJD(1) / 0.2e1 - t14 * mrSges(5,2) - t21 * mrSges(4,3) - t4 * mrSges(6,1) - t92 + (Ifges(6,6) + Ifges(4,6)) * t96 - (Ifges(5,3) + Ifges(6,2)) * t75 / 0.2e1;
t46 = t75 * t80;
t39 = -qJD(3) * mrSges(6,1) - mrSges(6,3) * t76;
t37 = t67 * qJD(1);
t36 = (mrSges(6,1) * t61 + mrSges(6,2) * t60) * qJD(1);
t35 = (-t61 * mrSges(5,1) - t60 * mrSges(5,3)) * qJD(1);
t31 = t81 * t60;
t22 = t65 * qJD(1);
t16 = t23 * qJD(1);
t13 = -qJD(3) * t32 - qJD(5) * t60;
t12 = t79 * t81 - t77;
t5 = t11 * qJD(1);
t2 = -qJD(5) * t76 + (t52 + (t38 - t74) * t61) * qJD(3);
t1 = (qJ(5) * t79 - t77) * qJD(1) + t7;
t24 = [t11 * t36 + t12 * t43 + t13 * t39 + t19 * t46 + t23 * t35 + m(6) * (-t1 * t32 + t11 * t4 + t12 * t6 + t13 * t3 + t19 * t5 - t2 * t31) + m(5) * (t15 * t23 + t16 * t30) + (-t16 * mrSges(5,1) + t5 * mrSges(6,1) + t7 * mrSges(5,2) + t17 * mrSges(4,3) - t1 * mrSges(6,3) + (t17 * t90 + t7 * t89) * t91 + (t71 * qJD(3) + (-m(4) * t20 + m(5) * t10 + t85) * t47 + (mrSges(4,2) * t72 - t30 * mrSges(5,3) + t31 * mrSges(6,3) - t61 * t70) * qJD(1) + t64) * qJD(3)) * t61 + (t5 * mrSges(6,2) - t16 * mrSges(5,3) - t2 * mrSges(6,3) + ((t89 + t90) * t91 + t101) * t18 + (t69 * qJD(3) + (-m(4) * t21 - m(5) * t14 - t83) * t47 + (t70 * t60 - t32 * mrSges(6,3) + t72 * mrSges(4,1) + t30 * mrSges(5,1) - t19 * mrSges(6,1) + (0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(4,1) + 0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(5,3) - 0.3e1 / 0.2e1 * Ifges(4,2)) * t61) * qJD(1) + t63) * qJD(3)) * t60; m(4) * (t17 * t60 - t87) + m(5) * (t60 * t7 - t87) + m(6) * (t1 * t60 - t2 * t61) + ((t43 + t83) * t61 + (t39 + t85) * t60 + m(4) * (-t20 * t60 + t21 * t61) + m(5) * (t10 * t60 + t14 * t61) + t98 + (t60 ^ 2 + t61 ^ 2) * qJD(1) * (-mrSges(4,3) - t86)) * qJD(3); -t2 * mrSges(6,1) - t17 * mrSges(4,2) + t1 * mrSges(6,2) + t7 * mrSges(5,3) - t22 * t36 - t37 * t35 - t9 * t39 - t8 * t43 - t85 * t21 - t83 * t20 + t102 * t18 + t84 * qJD(4) + ((Ifges(4,4) * t99 + (-t86 * qJ(4) + t69) * qJD(3) - t63 + t92) * t60 + (-t51 / 0.2e1 + (Ifges(6,4) / 0.2e1 + Ifges(5,5) / 0.2e1) * t75 + (-pkin(3) * mrSges(5,2) + mrSges(6,3) * t88 + t71) * qJD(3) + (-Ifges(6,1) / 0.2e1 + Ifges(4,2) / 0.2e1 - Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 - Ifges(5,1) / 0.2e1) * t76 - t64) * t61) * qJD(1) + (t1 * qJ(4) - t2 * t88 - t4 * t22 - t3 * t9 + t94 * t6) * m(6) + (-t18 * pkin(3) + t7 * qJ(4) - t10 * t21 + t93 * t14 - t15 * t37) * m(5); -t84 * qJD(3) + ((t35 - t36) * t60 + t86 * t61 * qJD(3)) * qJD(1) + (-t6 * qJD(3) - t4 * t76 + t2) * m(6) + (-t14 * qJD(3) + t15 * t76 + t18) * m(5); t46 + m(6) * t5 + (-mrSges(6,1) * t79 + t60 * t39 + t61 * t43 + t98) * qJD(1);];
tauc = t24(:);
