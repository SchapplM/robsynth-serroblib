% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PRPRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d4,d5,theta1]';
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
% Datum: 2019-12-31 17:40
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PRPRR9_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPRR9_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:39:37
% EndTime: 2019-12-31 17:39:39
% DurationCPUTime: 0.57s
% Computational Cost: add. (791->122), mult. (1377->175), div. (0->0), fcn. (555->4), ass. (0->62)
t35 = -qJD(2) + qJD(4);
t85 = t35 / 0.2e1;
t91 = -qJD(5) / 0.2e1;
t36 = sin(qJ(5));
t69 = t35 * t36;
t25 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t69;
t38 = cos(qJ(5));
t68 = t35 * t38;
t26 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t40 = -pkin(2) - pkin(3);
t30 = qJD(2) * t40 + qJD(3);
t37 = sin(qJ(4));
t39 = cos(qJ(4));
t59 = qJD(2) * qJ(3);
t17 = t30 * t37 + t39 * t59;
t8 = pkin(7) * t35 + t17;
t3 = -qJD(1) * t38 - t36 * t8;
t47 = qJD(1) * t36 - t38 * t8;
t51 = t3 * t36 + t38 * t47;
t90 = m(6) * t51 + t36 * t25 - t38 * t26;
t89 = (m(4) * qJ(3) + mrSges(4,3)) * qJD(2);
t88 = -Ifges(6,1) / 0.2e1;
t71 = Ifges(6,4) * t38;
t87 = (-Ifges(6,2) * t36 + t71) * t91;
t86 = -Ifges(6,4) * t68 / 0.2e1;
t72 = Ifges(6,4) * t36;
t84 = (Ifges(6,2) * t38 + t72) * t85;
t66 = t39 * qJ(3) + t37 * t40;
t60 = qJ(3) * qJD(4);
t61 = t39 * qJD(3);
t63 = qJD(4) * t30;
t5 = t39 * t63 + (-t37 * t60 + t61) * qJD(2);
t1 = qJD(5) * t3 + t38 * t5;
t2 = qJD(5) * t47 - t36 * t5;
t82 = t1 * t38 - t2 * t36;
t19 = (mrSges(6,1) * t36 + mrSges(6,2) * t38) * qJD(5);
t27 = -mrSges(6,1) * t38 + mrSges(6,2) * t36;
t62 = t37 * qJD(3);
t6 = t37 * t63 + (t39 * t60 + t62) * qJD(2);
t16 = t30 * t39 - t37 * t59;
t7 = -pkin(4) * t35 - t16;
t81 = -t5 * mrSges(5,2) + t7 * t19 + (t27 - mrSges(5,1)) * t6;
t42 = m(6) * (-t3 * t38 + t36 * t47) - t38 * t25 - t36 * t26;
t65 = Ifges(6,5) * qJD(5);
t53 = t3 * mrSges(6,3) + t69 * t88 + t86 - t65 / 0.2e1;
t64 = Ifges(6,6) * qJD(5);
t56 = -t47 * mrSges(6,3) + t64 / 0.2e1 + t84;
t80 = (Ifges(6,5) * t38 - Ifges(6,6) * t36) * t91 - ((Ifges(6,1) * t36 + t71) * t85 - t53) * t38 + (t84 + t56) * t36;
t79 = -m(5) * t17 + t35 * mrSges(5,2) + t90;
t18 = t27 * t35;
t57 = t35 * mrSges(5,1) - t18;
t55 = -t1 * mrSges(6,3) + t35 * t87;
t22 = (Ifges(6,1) * t38 - t72) * qJD(5);
t54 = -t2 * mrSges(6,3) + t22 * t85;
t48 = -qJ(3) * t37 + t39 * t40;
t46 = -t36 * t22 / 0.2e1 + t38 * t87;
t24 = -pkin(7) + t66;
t23 = pkin(4) - t48;
t15 = qJD(4) * t66 + t62;
t14 = qJD(4) * t48 + t61;
t11 = t35 * t19;
t4 = [m(6) * (-t1 * t36 - t2 * t38) + ((t36 ^ 2 + t38 ^ 2) * t35 * mrSges(6,3) + t90) * qJD(5); t23 * t11 + t15 * t18 + 0.2e1 * qJD(3) * t89 + (t14 * t26 + t55) * t38 + (-t14 * t25 - t54) * t36 + m(5) * (t17 * t14 - t16 * t15 - t48 * t6 + t5 * t66) + m(6) * (-t14 * t51 + t7 * t15 + t6 * t23 + t24 * t82) + (-t15 * mrSges(5,1) - t14 * mrSges(5,2) + t46) * t35 + (t42 * t24 + t80) * qJD(5) - t81; (-t11 + (-m(5) - m(6)) * t6 - t79 * qJD(4)) * t39 + (-t57 * qJD(4) + m(5) * (-qJD(4) * t16 + t5) + m(6) * (qJD(4) * t7 + t82) + t42 * qJD(5)) * t37 + ((m(5) * t16 - m(6) * t7 + t57) * t37 + t79 * t39 - t89) * qJD(2); -pkin(4) * t11 - t17 * t18 + (-t16 * t26 - t55) * t38 + (t16 * t25 + t54) * t36 + (t17 * mrSges(5,1) + t16 * mrSges(5,2) - t46) * t35 + (t42 * pkin(7) - t80) * qJD(5) + (-t6 * pkin(4) + pkin(7) * t82 + t51 * t16 - t7 * t17) * m(6) + t81; t2 * mrSges(6,1) - t1 * mrSges(6,2) - t47 * t25 - t3 * t26 + ((t65 / 0.2e1 - t7 * mrSges(6,2) + t86 + t53) * t38 + (-t64 / 0.2e1 - t7 * mrSges(6,1) + (t72 / 0.2e1 + (t88 + Ifges(6,2) / 0.2e1) * t38) * t35 + t56) * t36) * t35;];
tauc = t4(:);
