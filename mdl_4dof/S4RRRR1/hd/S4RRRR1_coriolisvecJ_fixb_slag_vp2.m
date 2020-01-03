% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,d3,d4]';
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
% Datum: 2019-12-31 17:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RRRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:22:08
% EndTime: 2019-12-31 17:22:09
% DurationCPUTime: 0.56s
% Computational Cost: add. (869->110), mult. (1738->176), div. (0->0), fcn. (820->6), ass. (0->70)
t45 = sin(qJ(4));
t48 = cos(qJ(4));
t102 = mrSges(5,1) * t48 - mrSges(5,2) * t45 + mrSges(4,1);
t47 = sin(qJ(2));
t50 = cos(qJ(2));
t101 = mrSges(3,1) * t47 + mrSges(3,2) * t50;
t42 = qJD(1) + qJD(2);
t82 = pkin(1) * qJD(1);
t34 = pkin(2) * t42 + t50 * t82;
t46 = sin(qJ(3));
t49 = cos(qJ(3));
t74 = t47 * t82;
t18 = t34 * t46 + t49 * t74;
t41 = qJD(3) + t42;
t12 = pkin(7) * t41 + t18;
t85 = t46 * t47;
t63 = t49 * t50 - t85;
t81 = qJD(3) * t46;
t53 = (qJD(2) * t63 - t47 * t81) * pkin(1);
t80 = qJD(3) * t49;
t6 = qJD(1) * t53 + t34 * t80;
t77 = qJD(4) * t48;
t3 = -t12 * t77 - t45 * t6;
t78 = qJD(4) * t45;
t2 = -t12 * t78 + t48 * t6;
t98 = t2 * t48;
t68 = -t3 * t45 + t98;
t73 = t12 * (t45 ^ 2 + t48 ^ 2);
t100 = t45 / 0.2e1;
t94 = Ifges(5,1) * t45;
t93 = Ifges(5,4) * t45;
t91 = Ifges(5,5) * t48;
t90 = Ifges(5,2) * t48;
t89 = Ifges(5,6) * t45;
t87 = t41 * t48;
t86 = t45 * mrSges(5,3);
t84 = t47 * t49;
t40 = pkin(1) * t50 + pkin(2);
t83 = pkin(1) * t84 + t46 * t40;
t28 = pkin(7) + t83;
t79 = qJD(4) * t28;
t71 = t77 / 0.2e1;
t70 = t102 * t41;
t66 = mrSges(5,1) * t45 + mrSges(5,2) * t48;
t32 = qJD(4) * mrSges(5,1) - t41 * t86;
t33 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t87;
t65 = t32 * t48 + t33 * t45;
t64 = t46 * t50 + t84;
t62 = -pkin(1) * t85 + t40 * t49;
t61 = t45 * (Ifges(5,1) * t48 - t93);
t60 = (t90 + t93) * t41;
t59 = t66 * qJD(4);
t58 = t65 * qJD(4);
t17 = t34 * t49 - t46 * t74;
t57 = t41 * mrSges(4,2) + t45 * t32 - t48 * t33;
t11 = -pkin(3) * t41 - t17;
t19 = Ifges(5,6) * qJD(4) + t60;
t35 = Ifges(5,4) * t87;
t20 = Ifges(5,5) * qJD(4) + t41 * t94 + t35;
t52 = (qJD(2) * t64 + t47 * t80) * pkin(1);
t7 = qJD(1) * t52 + t34 * t81;
t54 = mrSges(5,3) * t98 + t11 * t59 + t20 * t71 + qJD(4) ^ 2 * (-t89 + t91) / 0.2e1 - (t19 + t60) * t78 / 0.2e1 - t102 * t7 + (qJD(4) * t61 + (0.3e1 * Ifges(5,4) * t48 - 0.2e1 * Ifges(5,2) * t45 + t94) * t71) * t41;
t51 = -t6 * mrSges(4,2) - t3 * t86 + t54;
t30 = t63 * t82;
t29 = t64 * t82;
t27 = -pkin(3) - t62;
t21 = t41 * t59;
t10 = t40 * t81 + t52;
t9 = t40 * t80 + t53;
t1 = [t54 + m(5) * (t10 * t11 + t27 * t7 + t28 * t68 + t73 * t9) + (-t32 * t79 + t33 * t9) * t48 + (-t3 * mrSges(5,3) - t32 * t9 - t33 * t79) * t45 - t70 * t10 + m(4) * (-t17 * t10 + t18 * t9 + t6 * t83 - t62 * t7) + (-t41 * t9 - t6) * mrSges(4,2) + t27 * t21 + t101 * pkin(1) * qJD(2) * (-qJD(1) - t42); t51 + (m(4) * (t46 * t6 - t49 * t7) + (-t70 * t46 - t57 * t49 + m(4) * (-t17 * t46 + t18 * t49) + m(5) * (t11 * t46 + t49 * t73)) * qJD(3)) * pkin(2) - m(5) * (t11 * t29 + t30 * t73) + t70 * t29 + t57 * t30 - m(4) * (-t17 * t29 + t18 * t30) + t101 * t82 * (-qJD(2) + t42) + (m(5) * t7 + t21) * (-pkin(2) * t49 - pkin(3)) + (m(5) * t68 - t58) * (pkin(2) * t46 + pkin(7)); -pkin(3) * t21 - pkin(7) * t58 + t57 * t17 + t70 * t18 + t51 + (-t7 * pkin(3) + pkin(7) * t68 - t11 * t18 - t17 * t73) * m(5); t3 * mrSges(5,1) - t2 * mrSges(5,2) + t65 * t12 + (-t11 * t66 + t19 * t100 + (-t61 / 0.2e1 + t90 * t100) * t41 + (t91 / 0.2e1 - t89 / 0.2e1) * qJD(4) - (t20 + t35) * t48 / 0.2e1) * t41;];
tauc = t1(:);
