% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d4,theta3]';
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
% Datum: 2019-12-31 17:01
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRPR3_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPR3_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:01:28
% EndTime: 2019-12-31 17:01:29
% DurationCPUTime: 0.51s
% Computational Cost: add. (559->110), mult. (1257->176), div. (0->0), fcn. (622->6), ass. (0->63)
t47 = sin(qJ(4));
t49 = cos(qJ(4));
t89 = t49 * mrSges(5,1) - t47 * mrSges(5,2) + mrSges(4,1);
t48 = sin(qJ(2));
t50 = cos(qJ(2));
t88 = mrSges(3,1) * t48 + mrSges(3,2) * t50;
t46 = cos(pkin(7));
t45 = sin(pkin(7));
t80 = t45 * t48;
t56 = pkin(1) * (t46 * t50 - t80);
t31 = qJD(2) * t56;
t18 = qJD(1) * t31;
t44 = qJD(1) + qJD(2);
t76 = pkin(1) * qJD(1);
t37 = t44 * pkin(2) + t50 * t76;
t67 = t48 * t76;
t13 = t45 * t37 + t46 * t67;
t9 = t44 * pkin(6) + t13;
t5 = t49 * qJD(3) - t47 * t9;
t75 = qJD(4) * t5;
t2 = t49 * t18 + t75;
t87 = t2 * t49;
t84 = Ifges(5,1) * t47;
t83 = Ifges(5,4) * t47;
t81 = t44 * t49;
t79 = t46 * t48;
t78 = t47 * mrSges(5,3);
t43 = t50 * pkin(1) + pkin(2);
t77 = pkin(1) * t79 + t45 * t43;
t6 = t47 * qJD(3) + t49 * t9;
t74 = qJD(4) * t6;
t73 = Ifges(5,5) * qJD(4);
t72 = Ifges(5,6) * qJD(4);
t27 = pkin(6) + t77;
t71 = qJD(4) * t27;
t70 = qJD(4) * t47;
t65 = qJD(4) * t49 / 0.2e1;
t64 = t89 * t44;
t62 = -t6 * t47 - t5 * t49;
t61 = -t47 * t5 + t49 * t6;
t35 = qJD(4) * mrSges(5,1) - t44 * t78;
t36 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t81;
t59 = -t47 * t35 + t49 * t36;
t58 = -pkin(1) * t80 + t46 * t43;
t57 = pkin(1) * (t45 * t50 + t79);
t55 = (Ifges(5,2) * t49 + t83) * t44;
t54 = (mrSges(5,1) * t47 + mrSges(5,2) * t49) * qJD(4);
t12 = t46 * t37 - t45 * t67;
t29 = qJD(2) * t57;
t17 = qJD(1) * t29;
t21 = t55 + t72;
t39 = Ifges(5,4) * t81;
t22 = t44 * t84 + t39 + t73;
t8 = -t44 * pkin(3) - t12;
t51 = mrSges(5,3) * t87 + t22 * t65 + t8 * t54 + qJD(4) ^ 2 * (Ifges(5,5) * t49 - Ifges(5,6) * t47) / 0.2e1 - (t21 + t55) * t70 / 0.2e1 - t89 * t17 + ((Ifges(5,1) * t49 - t83) * t70 + (0.3e1 * Ifges(5,4) * t49 - 0.2e1 * Ifges(5,2) * t47 + t84) * t65) * t44;
t42 = -t46 * pkin(2) - pkin(3);
t41 = t45 * pkin(2) + pkin(6);
t30 = qJD(1) * t56;
t28 = qJD(1) * t57;
t26 = -pkin(3) - t58;
t23 = t44 * t54;
t3 = -t47 * t18 - t74;
t1 = [t51 + (-mrSges(5,3) * t75 - t35 * t71 + t31 * t36 + m(5) * (t2 * t27 + t31 * t6 - t5 * t71)) * t49 + (-t36 * t71 - t31 * t35 + m(5) * (-t27 * t3 - t31 * t5 - t6 * t71) + (-t3 - t74) * mrSges(5,3)) * t47 - t64 * t29 + m(5) * (t17 * t26 + t8 * t29) + m(4) * (-t12 * t29 + t13 * t31 - t17 * t58 + t18 * t77) + (-t31 * t44 - t18) * mrSges(4,2) + t26 * t23 + t88 * pkin(1) * qJD(2) * (-qJD(1) - t44); t51 + ((-t35 * t49 - t36 * t47) * t41 + t62 * mrSges(5,3)) * qJD(4) + (t44 * mrSges(4,2) - t59) * t30 + t64 * t28 + t42 * t23 - t18 * mrSges(4,2) - t3 * t78 + t88 * t76 * (-qJD(2) + t44) + ((-t17 * t46 + t18 * t45) * pkin(2) + t12 * t28 - t13 * t30) * m(4) + (-t8 * t28 - t30 * t61 + (qJD(4) * t62 - t3 * t47 + t87) * t41 + t17 * t42) * m(5); m(5) * (t2 * t47 + t3 * t49) + (m(5) * t61 + (-t47 ^ 2 - t49 ^ 2) * t44 * mrSges(5,3) + t59) * qJD(4); t3 * mrSges(5,1) - t2 * mrSges(5,2) + t6 * t35 - t5 * t36 + ((t73 / 0.2e1 - t8 * mrSges(5,2) - t22 / 0.2e1 - t39 / 0.2e1 + t5 * mrSges(5,3)) * t49 + (-t72 / 0.2e1 - t8 * mrSges(5,1) + t21 / 0.2e1 + t6 * mrSges(5,3) + (t83 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t49) * t44) * t47) * t44;];
tauc = t1(:);
