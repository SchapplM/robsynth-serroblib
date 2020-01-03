% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:39
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:39:43
% EndTime: 2019-12-31 16:39:44
% DurationCPUTime: 0.38s
% Computational Cost: add. (330->84), mult. (725->134), div. (0->0), fcn. (311->4), ass. (0->52)
t24 = sin(qJ(4));
t25 = cos(qJ(4));
t26 = -pkin(1) - pkin(2);
t18 = t26 * qJD(1) + qJD(2);
t22 = sin(pkin(6));
t23 = cos(pkin(6));
t50 = qJD(1) * qJ(2);
t8 = t22 * t18 + t23 * t50;
t6 = -qJD(1) * pkin(5) + t8;
t3 = t25 * qJD(3) - t24 * t6;
t4 = t24 * qJD(3) + t25 * t6;
t36 = -t24 * t3 + t25 * t4;
t54 = qJD(1) * mrSges(5,3);
t16 = qJD(4) * mrSges(5,1) + t24 * t54;
t17 = -qJD(4) * mrSges(5,2) - t25 * t54;
t70 = -t24 * t16 + t25 * t17;
t76 = -m(5) * t36 - t70;
t49 = qJD(1) * qJD(2);
t42 = t23 * t49;
t1 = qJD(4) * t3 + t25 * t42;
t2 = -qJD(4) * t4 - t24 * t42;
t39 = t1 * t25 - t2 * t24;
t51 = qJD(4) * t25;
t52 = qJD(4) * t24;
t75 = t3 * t51 + t4 * t52 - t39;
t73 = -m(3) * qJ(2) - mrSges(3,3);
t72 = -t16 * t51 - t17 * t52;
t64 = Ifges(5,4) * t24;
t63 = Ifges(5,4) * t25;
t62 = Ifges(5,5) * t25;
t61 = Ifges(5,6) * t24;
t60 = t22 * qJD(1) * (t25 * mrSges(5,1) - t24 * mrSges(5,2));
t59 = t23 * t18;
t56 = t23 * qJ(2) + t22 * t26;
t55 = t22 * qJ(2);
t53 = qJD(2) * t23;
t38 = -(-t22 * t50 + t59) * t22 + t8 * t23;
t37 = -t24 * t4 - t25 * t3;
t35 = -mrSges(5,1) * t24 - mrSges(5,2) * t25;
t33 = t23 * t26 - t55;
t32 = t24 * (-Ifges(5,1) * t25 + t64);
t31 = t25 * (Ifges(5,2) * t24 - t63);
t30 = t35 * qJD(4);
t29 = (-Ifges(5,1) * t24 - t63) * qJD(1);
t28 = (-Ifges(5,2) * t25 - t64) * qJD(1);
t27 = qJD(1) ^ 2;
t12 = pkin(3) - t33;
t11 = qJD(1) * t30;
t10 = Ifges(5,5) * qJD(4) + t29;
t9 = Ifges(5,6) * qJD(4) + t28;
t5 = -t59 + (pkin(3) + t55) * qJD(1);
t7 = [qJD(4) ^ 2 * (t61 - t62) / 0.2e1 + t12 * t11 + 0.2e1 * mrSges(4,2) * t42 + t5 * t30 + t70 * t53 + (-t32 - t31) * qJD(4) * qJD(1) + (t28 + t9) * t52 / 0.2e1 - (t29 + t10) * t51 / 0.2e1 + (m(5) * (t37 * qJD(4) + t39) + t72) * (-pkin(5) + t56) + (m(5) * (t36 * t23 + (qJD(1) * t12 + t5) * t22) + m(4) * ((-t22 * t33 + t23 * t56) * qJD(1) + t38) + 0.2e1 * t60) * qJD(2) + 0.2e1 * (t22 * mrSges(4,1) - t73) * t49 + t75 * mrSges(5,3); t73 * t27 + (-m(4) * t38 - t60) * qJD(1) + (-t27 * mrSges(4,1) + ((-t5 - t53) * qJD(1) - t75) * m(5) + t72) * t22 + (-t27 * mrSges(4,2) + t76 * qJD(1) - t11) * t23; m(5) * (t1 * t24 + t2 * t25) + ((t24 ^ 2 + t25 ^ 2) * t54 - t76) * qJD(4); t2 * mrSges(5,1) - t1 * mrSges(5,2) + t4 * t16 - t3 * t17 + (-t5 * t35 + t25 * t10 / 0.2e1 - t24 * t9 / 0.2e1 + (t32 / 0.2e1 + t31 / 0.2e1) * qJD(1) + t37 * mrSges(5,3) + (-t62 / 0.2e1 + t61 / 0.2e1) * qJD(4)) * qJD(1);];
tauc = t7(:);
