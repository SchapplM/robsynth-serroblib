% Calculate vector of centrifugal and Coriolis load on the joints for
% S5RPPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2]';
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
% Datum: 2019-12-31 17:56
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5RPPRR5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:56:29
% EndTime: 2019-12-31 17:56:31
% DurationCPUTime: 0.62s
% Computational Cost: add. (1023->124), mult. (1802->178), div. (0->0), fcn. (786->6), ass. (0->62)
t38 = -qJD(1) + qJD(4);
t90 = t38 / 0.2e1;
t41 = sin(qJ(5));
t69 = t38 * t41;
t26 = qJD(5) * mrSges(6,1) - mrSges(6,3) * t69;
t43 = cos(qJ(5));
t68 = t38 * t43;
t27 = -qJD(5) * mrSges(6,2) + mrSges(6,3) * t68;
t35 = -cos(pkin(8)) * pkin(1) - pkin(2) - pkin(3);
t21 = t35 * qJD(1) + qJD(3);
t36 = sin(pkin(8)) * pkin(1) + qJ(3);
t29 = qJD(1) * t36;
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t10 = t42 * t21 + t44 * t29;
t8 = t38 * pkin(7) + t10;
t3 = -t43 * qJD(2) - t41 * t8;
t50 = t41 * qJD(2) - t43 * t8;
t54 = t3 * t41 + t43 * t50;
t89 = m(6) * t54 + t41 * t26 - t43 * t27;
t88 = -Ifges(6,1) / 0.2e1;
t87 = -Ifges(6,4) * t68 / 0.2e1;
t86 = -t38 / 0.2e1;
t72 = Ifges(6,4) * t41;
t85 = (Ifges(6,2) * t43 + t72) * t90;
t67 = t42 * t35 + t44 * t36;
t61 = qJD(1) * qJD(3);
t9 = t44 * t21 - t42 * t29;
t65 = qJD(4) * t9;
t5 = t44 * t61 + t65;
t1 = t3 * qJD(5) + t43 * t5;
t2 = t50 * qJD(5) - t41 * t5;
t83 = t1 * t43 - t2 * t41;
t22 = (mrSges(6,1) * t41 + mrSges(6,2) * t43) * qJD(5);
t31 = -t43 * mrSges(6,1) + t41 * mrSges(6,2);
t62 = qJD(4) * t10;
t6 = t42 * t61 + t62;
t7 = -t38 * pkin(4) - t9;
t82 = -t5 * mrSges(5,2) + t7 * t22 + (t31 - mrSges(5,1)) * t6;
t81 = -t38 * mrSges(5,2) - t89;
t45 = m(6) * (-t3 * t43 + t41 * t50) - t43 * t26 - t41 * t27;
t63 = Ifges(6,6) * qJD(5);
t56 = -t50 * mrSges(6,3) + t63 / 0.2e1 + t85;
t64 = Ifges(6,5) * qJD(5);
t57 = t3 * mrSges(6,3) + t69 * t88 + t87 - t64 / 0.2e1;
t71 = Ifges(6,4) * t43;
t80 = -qJD(5) * (Ifges(6,5) * t43 - Ifges(6,6) * t41) / 0.2e1 - ((Ifges(6,1) * t41 + t71) * t90 - t57) * t43 + (t85 + t56) * t41;
t79 = m(4) * t29;
t20 = t31 * t38;
t60 = t38 * mrSges(5,1) - t20;
t24 = (-Ifges(6,2) * t41 + t71) * qJD(5);
t59 = -t1 * mrSges(6,3) + t24 * t86;
t25 = (Ifges(6,1) * t43 - t72) * qJD(5);
t58 = t2 * mrSges(6,3) + t25 * t86;
t51 = t44 * t35 - t42 * t36;
t49 = -t41 * t25 / 0.2e1 - t43 * t24 / 0.2e1;
t17 = t38 * t22;
t14 = -pkin(7) + t67;
t13 = pkin(4) - t51;
t12 = t42 * qJD(3) + t67 * qJD(4);
t11 = t44 * qJD(3) + t51 * qJD(4);
t4 = [0.2e1 * mrSges(4,3) * t61 + t12 * t20 + t13 * t17 + (t11 * t27 + t59) * t43 + (-t11 * t26 + t58) * t41 + m(6) * (-t54 * t11 + t7 * t12 + t6 * t13 + t83 * t14) + m(5) * (t10 * t11 - t9 * t12 + t5 * t67 - t6 * t51) + 0.2e1 * qJD(3) * t79 + (-t12 * mrSges(5,1) - t11 * mrSges(5,2) + t49) * t38 + (t45 * t14 + t80) * qJD(5) - t82; m(6) * (-t1 * t41 - t2 * t43) + ((t41 ^ 2 + t43 ^ 2) * t38 * mrSges(6,3) + t89) * qJD(5); (-mrSges(4,3) * qJD(1) - t79) * qJD(1) + (-t60 * qJD(4) + m(5) * (t5 - t65) + m(6) * (qJD(4) * t7 + t83) + (m(5) * t9 - m(6) * t7 + t60) * qJD(1) + t45 * qJD(5)) * t42 + (-t17 + m(5) * (-t6 + t62) - m(6) * t6 + (-m(5) * t10 - t81) * qJD(1) + t81 * qJD(4)) * t44; -pkin(4) * t17 - t10 * t20 + (-t9 * t27 - t59) * t43 + (t9 * t26 - t58) * t41 + (t10 * mrSges(5,1) + t9 * mrSges(5,2) - t49) * t38 + (t45 * pkin(7) - t80) * qJD(5) + (-t6 * pkin(4) + t83 * pkin(7) - t7 * t10 + t54 * t9) * m(6) + t82; t2 * mrSges(6,1) - t1 * mrSges(6,2) - t50 * t26 - t3 * t27 + ((t64 / 0.2e1 - t7 * mrSges(6,2) + t87 + t57) * t43 + (-t63 / 0.2e1 - t7 * mrSges(6,1) + (t72 / 0.2e1 + (t88 + Ifges(6,2) / 0.2e1) * t43) * t38 + t56) * t41) * t38;];
tauc = t4(:);
