% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPRPR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
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
% tauc [5x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 18:18
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPRPR6_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR6_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:17:42
% EndTime: 2019-12-31 18:17:43
% DurationCPUTime: 0.56s
% Computational Cost: add. (830->115), mult. (1622->163), div. (0->0), fcn. (750->6), ass. (0->67)
t48 = cos(qJ(5));
t97 = -t48 / 0.2e1;
t47 = sin(qJ(3));
t49 = cos(qJ(3));
t42 = cos(pkin(8)) * pkin(1) + pkin(2);
t58 = t42 * qJD(1);
t94 = pkin(1) * sin(pkin(8));
t77 = qJD(1) * t94;
t20 = t47 * t77 - t49 * t58;
t57 = t20 + qJD(4);
t46 = sin(qJ(5));
t96 = (-Ifges(6,5) * t46 / 0.2e1 + Ifges(6,6) * t97) * qJD(5);
t83 = mrSges(4,1) - mrSges(5,2);
t21 = t47 * t58 + t49 * t77;
t44 = qJD(1) + qJD(3);
t11 = t44 * qJ(4) + t21;
t55 = qJD(3) * t58;
t76 = qJD(3) * t94;
t72 = qJD(1) * t76;
t16 = -t47 * t72 + t49 * t55;
t9 = -qJD(4) * t44 - t16;
t95 = -t9 * qJ(4) + t11 * t57;
t50 = -pkin(3) - pkin(7);
t8 = t44 * t50 + t57;
t4 = -qJD(2) * t46 + t48 * t8;
t5 = qJD(2) * t48 + t46 * t8;
t70 = t4 * t48 + t46 * t5;
t17 = t47 * t55 + t49 * t72;
t1 = qJD(5) * t4 + t46 * t17;
t93 = t1 * mrSges(6,3);
t90 = Ifges(6,4) * t46;
t89 = Ifges(6,4) * t48;
t85 = t44 * mrSges(6,3);
t84 = t49 * t42;
t82 = qJD(5) * t5;
t73 = -t47 * t94 + t84;
t71 = -pkin(3) - t73;
t30 = -pkin(7) + t71;
t81 = qJD(5) * t30;
t79 = qJD(5) * t46;
t78 = qJD(5) * t50;
t69 = -t4 * t46 + t48 * t5;
t25 = qJD(3) * t84 - t47 * t76;
t22 = -qJD(4) - t25;
t63 = t42 * t47 + t49 * t94;
t33 = qJ(4) + t63;
t68 = -t11 * t22 - t9 * t33;
t2 = t48 * t17 - t82;
t67 = (-t2 - t82) * mrSges(6,3);
t66 = mrSges(6,1) * t48 - mrSges(6,2) * t46;
t65 = mrSges(6,1) * t46 + mrSges(6,2) * t48;
t62 = t46 * (-Ifges(6,2) * t48 - t90);
t61 = t48 * (-Ifges(6,1) * t46 - t89);
t60 = (Ifges(6,1) * t48 - t90) * t44;
t59 = (-Ifges(6,2) * t46 + t89) * t44;
t56 = t66 * qJD(5);
t52 = m(6) * (qJD(5) * t69 + t1 * t46 + t2 * t48);
t27 = Ifges(6,6) * qJD(5) + t59;
t28 = Ifges(6,5) * qJD(5) + t60;
t51 = t4 * mrSges(6,3) * t79 - t16 * mrSges(4,2) + t11 * t56 + (-mrSges(5,3) - t65) * t9 - (t28 + t60) * t79 / 0.2e1 - t83 * t17 + ((t61 - t62) * t44 + (t27 + t59) * t97 + t96) * qJD(5);
t39 = qJD(5) * mrSges(6,1) - t48 * t85;
t38 = -qJD(5) * mrSges(6,2) - t46 * t85;
t34 = t65 * t44;
t29 = t44 * t56;
t26 = t63 * qJD(3);
t10 = -pkin(3) * t44 + t57;
t3 = [t51 + (-t25 * mrSges(4,2) - t22 * mrSges(5,3) - t26 * t83) * t44 + m(6) * (t26 * t70 + t68) + t30 * t52 + (t26 * t39 + t38 * t81 + t67) * t48 + (t26 * t38 - t39 * t81 - t93) * t46 + m(4) * (t16 * t63 - t17 * t73 + t20 * t26 + t21 * t25) + m(5) * (t10 * t26 + t17 * t71 + t68) + t33 * t29 - t22 * t34; m(6) * (t1 * t48 - t2 * t46) + (-t46 * t38 - t48 * t39 - m(6) * t70 + (-t46 ^ 2 - t48 ^ 2) * t85) * qJD(5); t51 + t50 * t52 + (-t20 * mrSges(4,2) + mrSges(5,3) * t57 + t21 * t83) * t44 + (-t21 * t39 + t38 * t78 + t67) * t48 + (-t21 * t38 - t39 * t78 - t93) * t46 + t57 * t34 + qJ(4) * t29 + (-t21 * t70 + t95) * m(6) + (-t17 * pkin(3) - t10 * t21 + t95) * m(5); (t48 * t38 - t46 * t39) * qJD(5) + m(5) * t17 + t52 + (-mrSges(5,3) * t44 - t34 + (-m(5) - m(6)) * t11) * t44; t2 * mrSges(6,1) - t1 * mrSges(6,2) - t4 * t38 + t5 * t39 + (-t11 * t66 + t46 * t28 / 0.2e1 + t48 * t27 / 0.2e1 + (-t61 / 0.2e1 + t62 / 0.2e1) * t44 + t69 * mrSges(6,3) + t96) * t44;];
tauc = t3(:);
