% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta3,theta4]';
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
% Datum: 2019-12-31 17:46
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPPR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:46:21
% EndTime: 2019-12-31 17:46:23
% DurationCPUTime: 0.83s
% Computational Cost: add. (880->141), mult. (1870->218), div. (0->0), fcn. (1105->6), ass. (0->77)
t64 = cos(pkin(7));
t82 = t64 * qJD(2);
t50 = -qJD(4) + t82;
t47 = t50 * qJD(1);
t61 = sin(pkin(8));
t63 = cos(pkin(8));
t86 = t61 ^ 2 + t63 ^ 2;
t77 = t86 * t47;
t66 = sin(qJ(5));
t67 = cos(qJ(5));
t42 = t61 * t66 - t67 * t63;
t109 = t42 * qJD(5);
t108 = t109 / 0.2e1;
t43 = t61 * t67 + t63 * t66;
t41 = t43 * qJD(5);
t107 = t41 / 0.2e1;
t106 = mrSges(5,3) * t86;
t32 = qJD(1) * t109;
t105 = t42 * t32;
t33 = qJD(1) * t41;
t104 = t43 * t33;
t102 = -m(3) * qJ(2) - mrSges(3,3);
t62 = sin(pkin(7));
t84 = qJD(1) * t64;
t101 = t109 * t62 + t43 * t84;
t30 = t43 * t62;
t100 = -qJD(5) * t30 + t42 * t84;
t68 = -pkin(1) - pkin(2);
t48 = qJD(1) * t68 + qJD(2);
t81 = qJ(2) * qJD(1);
t35 = t62 * t48 + t64 * t81;
t20 = -qJD(1) * qJ(4) + t35;
t16 = t61 * qJD(3) + t63 * t20;
t56 = t63 * qJD(3);
t72 = t16 * t63 - (-t20 * t61 + t56) * t61;
t39 = t43 * qJD(1);
t97 = -t39 / 0.2e1;
t96 = m(5) * t62;
t87 = t64 * qJ(2) + t62 * t68;
t44 = -qJ(4) + t87;
t95 = pkin(6) - t44;
t93 = Ifges(6,4) * t39;
t85 = qJD(1) * t63;
t83 = t62 * qJD(2);
t80 = qJD(1) * qJD(2);
t75 = -t82 / 0.2e1;
t9 = -t33 * mrSges(6,1) + t32 * mrSges(6,2);
t34 = t48 * t64 - t62 * t81;
t74 = -t62 * qJ(2) + t64 * t68;
t73 = pkin(3) - t74;
t12 = t56 + (pkin(6) * qJD(1) - t20) * t61;
t13 = -pkin(6) * t85 + t16;
t3 = t12 * t67 - t13 * t66;
t4 = t12 * t66 + t13 * t67;
t25 = t95 * t61;
t26 = t95 * t63;
t7 = t25 * t67 + t26 * t66;
t8 = t25 * t66 - t26 * t67;
t71 = -t34 * t62 + t35 * t64;
t19 = qJD(1) * pkin(3) + qJD(4) - t34;
t45 = qJD(1) * (t63 * mrSges(5,1) - t61 * mrSges(5,2));
t69 = qJD(1) ^ 2;
t38 = t42 * qJD(1);
t37 = pkin(4) * t63 + t73;
t36 = Ifges(6,4) * t38;
t31 = t42 * t62;
t28 = qJD(5) * mrSges(6,1) + mrSges(6,3) * t39;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t38;
t17 = pkin(4) * t85 + t19;
t14 = -mrSges(6,1) * t38 - mrSges(6,2) * t39;
t11 = -Ifges(6,1) * t39 + Ifges(6,5) * qJD(5) + t36;
t10 = Ifges(6,2) * t38 + Ifges(6,6) * qJD(5) - t93;
t6 = -qJD(5) * t8 - t43 * t50;
t5 = qJD(5) * t7 - t42 * t50;
t2 = -qJD(5) * t4 - t43 * t47;
t1 = qJD(5) * t3 - t42 * t47;
t15 = [m(5) * (t44 * t77 + t72 * t50) + m(6) * (t1 * t8 + t2 * t7 + t3 * t6 + t4 * t5) + m(4) * ((-t62 * t74 + t64 * t87) * qJD(1) + t71) * qJD(2) + t11 * t108 + t17 * (-t41 * mrSges(6,1) + mrSges(6,2) * t109) + qJD(5) * (Ifges(6,5) * t109 + Ifges(6,6) * t41) / 0.2e1 + t10 * t107 + t37 * t9 + t5 * t27 + t6 * t28 + (m(5) * (qJD(1) * t73 + t19) + m(6) * (qJD(1) * t37 + t17) + 0.2e1 * t45 + t14) * t83 + (-t42 * mrSges(6,1) - t43 * mrSges(6,2) + (2 * mrSges(4,1))) * t62 * t80 + (t38 * t107 + t42 * t33) * Ifges(6,2) + (t109 * t97 - t32 * t43) * Ifges(6,1) + 0.2e1 * (t64 * mrSges(4,2) - t102) * t80 + (t1 * t42 - t109 * t3 + t2 * t43 - t32 * t7 + t33 * t8 + t4 * t41) * mrSges(6,3) + (t38 * t108 + t41 * t97 - t104 + t105) * Ifges(6,4) - 0.2e1 * mrSges(5,3) * t77; -t64 * t9 + t101 * t28 + t100 * t27 + (t30 * t32 - t31 * t33) * mrSges(6,3) + (-t62 * mrSges(4,1) + (-mrSges(4,2) + t106) * t64 + t102) * t69 + t77 * t96 + (-t1 * t31 + t100 * t4 + t101 * t3 - t2 * t30) * m(6) + (-m(4) * t71 + 0.2e1 * t75 * t96 - m(5) * t72 * t64 + (-t14 - t45 - m(5) * t19 + 0.2e1 * (-t17 / 0.2e1 + t75) * m(6)) * t62) * qJD(1); m(6) * (t1 * t43 - t109 * t4 - t2 * t42 - t3 * t41) - t109 * t27 - t41 * t28 + (t104 + t105) * mrSges(6,3); -m(6) * (t3 * t39 + t38 * t4) - t38 * t27 - t39 * t28 - t69 * t106 + (m(6) * t83 + (t72 + t83) * m(5)) * qJD(1) + t9; Ifges(6,5) * t32 + Ifges(6,6) * t33 - t1 * mrSges(6,2) + t2 * mrSges(6,1) - t17 * (-mrSges(6,1) * t39 + mrSges(6,2) * t38) + t39 * (Ifges(6,1) * t38 + t93) / 0.2e1 + t10 * t97 - qJD(5) * (Ifges(6,5) * t38 + Ifges(6,6) * t39) / 0.2e1 - t3 * t27 + t4 * t28 + (t3 * t38 - t39 * t4) * mrSges(6,3) - (Ifges(6,2) * t39 + t11 + t36) * t38 / 0.2e1;];
tauc = t15(:);
