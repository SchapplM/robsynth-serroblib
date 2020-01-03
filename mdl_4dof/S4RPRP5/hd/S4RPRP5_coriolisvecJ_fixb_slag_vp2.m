% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [4x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 16:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:44:46
% EndTime: 2019-12-31 16:44:48
% DurationCPUTime: 0.95s
% Computational Cost: add. (650->152), mult. (1907->197), div. (0->0), fcn. (1212->4), ass. (0->75)
t92 = Ifges(4,1) + Ifges(5,1);
t91 = Ifges(5,4) + Ifges(4,5);
t93 = -mrSges(5,1) - mrSges(4,1);
t69 = Ifges(5,5) - Ifges(4,4);
t90 = -Ifges(4,6) + Ifges(5,6);
t51 = cos(pkin(6));
t79 = cos(qJ(3));
t60 = t79 * t51;
t57 = qJD(1) * t60;
t50 = sin(pkin(6));
t52 = sin(qJ(3));
t73 = t52 * t50;
t33 = qJD(1) * t73 - t57;
t31 = Ifges(4,4) * t33;
t39 = t79 * t50 + t52 * t51;
t34 = t39 * qJD(1);
t77 = Ifges(5,5) * t33;
t89 = t91 * qJD(3) + t92 * t34 - t31 + t77;
t68 = pkin(5) + qJ(2);
t43 = t68 * t50;
t40 = qJD(1) * t43;
t61 = t79 * t40;
t44 = t68 * t51;
t41 = qJD(1) * t44;
t74 = t52 * t41;
t18 = -t61 - t74;
t88 = -t18 + qJD(4);
t87 = (m(3) * qJ(2) + mrSges(3,3)) * (t50 ^ 2 + t51 ^ 2);
t86 = -t33 / 0.2e1;
t85 = t33 / 0.2e1;
t83 = t34 / 0.2e1;
t19 = -t52 * t40 + t79 * t41;
t54 = t39 * qJD(2);
t4 = qJD(1) * t54 + qJD(3) * t19;
t56 = -t79 * t43 - t52 * t44;
t81 = t56 * t4;
t80 = t4 * t39;
t78 = Ifges(4,4) * t34;
t76 = t33 * mrSges(4,3);
t75 = t34 * mrSges(4,3);
t72 = -qJD(3) / 0.2e1;
t71 = qJD(3) / 0.2e1;
t70 = mrSges(5,2) + mrSges(4,3);
t25 = -mrSges(5,2) * t33 + qJD(3) * mrSges(5,3);
t67 = -qJD(3) * mrSges(4,2) + t25 - t76;
t66 = -t34 * mrSges(5,2) - t93 * qJD(3) - t75;
t65 = qJD(2) * t57 - qJD(3) * t61;
t63 = qJD(1) * qJD(2);
t62 = -pkin(2) * t51 - pkin(1);
t59 = t50 * t63;
t21 = -t52 * t43 + t79 * t44;
t55 = t60 - t73;
t42 = qJD(1) * t62 + qJD(2);
t36 = t39 * qJD(3);
t35 = t55 * qJD(3);
t30 = Ifges(5,5) * t34;
t29 = qJD(1) * t36;
t28 = qJD(1) * t35;
t27 = t29 * mrSges(5,1);
t26 = t28 * mrSges(4,2);
t17 = -pkin(3) * t55 - qJ(4) * t39 + t62;
t16 = mrSges(5,1) * t33 - mrSges(5,3) * t34;
t15 = pkin(3) * t34 + qJ(4) * t33;
t14 = qJD(3) * qJ(4) + t19;
t13 = -qJD(3) * pkin(3) + t88;
t10 = -Ifges(4,2) * t33 + Ifges(4,6) * qJD(3) + t78;
t9 = Ifges(5,6) * qJD(3) + Ifges(5,3) * t33 + t30;
t8 = qJD(3) * t21 + t54;
t7 = qJD(2) * t55 + qJD(3) * t56;
t6 = pkin(3) * t33 - qJ(4) * t34 + t42;
t5 = pkin(3) * t36 - qJ(4) * t35 - qJD(4) * t39;
t3 = (-qJD(3) * t41 - t59) * t52 + t65;
t2 = -t52 * t59 + (qJD(4) - t74) * qJD(3) + t65;
t1 = pkin(3) * t29 - qJ(4) * t28 - qJD(4) * t34;
t11 = [t17 * t27 + t1 * (-mrSges(5,1) * t55 - t39 * mrSges(5,3)) + (-mrSges(5,3) * t17 + t39 * t92 - t55 * t69 - t56 * t70) * t28 + (t62 * mrSges(4,1) + t69 * t39 - (Ifges(5,3) + Ifges(4,2)) * t55 - t70 * t21) * t29 + 0.2e1 * t87 * t63 - t66 * t8 + t67 * t7 + t62 * t26 + t5 * t16 + (t3 * t55 + t80) * mrSges(4,3) + (t2 * t55 + t80) * mrSges(5,2) + m(5) * (t1 * t17 + t13 * t8 + t14 * t7 + t2 * t21 + t5 * t6 - t81) + m(4) * (-t18 * t8 + t19 * t7 + t3 * t21 - t81) + (Ifges(5,3) * t85 - Ifges(4,2) * t86 + t42 * mrSges(4,1) + t6 * mrSges(5,1) + t9 / 0.2e1 - t10 / 0.2e1 - t19 * mrSges(4,3) - t14 * mrSges(5,2) + t69 * t83 + t90 * t71) * t36 + (t89 / 0.2e1 + t42 * mrSges(4,2) + t13 * mrSges(5,2) - t18 * mrSges(4,3) - t6 * mrSges(5,3) + Ifges(4,4) * t86 + Ifges(5,5) * t85 + t91 * t71 + t92 * t83) * t35; t29 * mrSges(4,1) - t28 * mrSges(5,3) + t26 + t27 + t66 * t34 + t67 * t33 - m(4) * (-t18 * t34 - t19 * t33) - t87 * qJD(1) ^ 2 + (-t13 * t34 + t14 * t33 + t1) * m(5); -t42 * (mrSges(4,1) * t34 - mrSges(4,2) * t33) + t10 * t83 + (Ifges(5,3) * t34 - t77) * t86 - t6 * (t34 * mrSges(5,1) + t33 * mrSges(5,3)) + qJD(4) * t25 - t15 * t16 - t3 * mrSges(4,2) + t2 * mrSges(5,3) + t93 * t4 + t90 * t29 + t91 * t28 + (t66 + t75) * t19 + (-t67 - t76) * t18 + (-pkin(3) * t28 - qJ(4) * t29 + t13 * t33 + t14 * t34) * mrSges(5,2) + (-t33 * t91 + t34 * t90) * t72 + (-t4 * pkin(3) + t2 * qJ(4) - t13 * t19 + t14 * t88 - t6 * t15) * m(5) + (-Ifges(4,2) * t34 - t31 + t89) * t85 - (-t33 * t92 + t30 - t78 + t9) * t34 / 0.2e1; t28 * mrSges(5,2) - qJD(3) * t25 + t34 * t16 + 0.2e1 * (t4 / 0.2e1 + t14 * t72 + t6 * t83) * m(5);];
tauc = t11(:);
