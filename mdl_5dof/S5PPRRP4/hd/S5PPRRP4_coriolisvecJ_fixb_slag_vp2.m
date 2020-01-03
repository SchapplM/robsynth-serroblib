% Calculate vector of centrifugal and Coriolis load on the joints for
% S5PPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,theta1]';
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
% Datum: 2019-12-31 17:34
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S5PPRRP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:34:27
% EndTime: 2019-12-31 17:34:29
% DurationCPUTime: 0.85s
% Computational Cost: add. (464->155), mult. (1243->205), div. (0->0), fcn. (583->4), ass. (0->67)
t90 = -qJD(4) / 0.2e1;
t89 = Ifges(5,4) + Ifges(6,4);
t40 = sin(qJ(4));
t42 = cos(qJ(4));
t41 = sin(qJ(3));
t67 = qJD(2) * t41;
t29 = qJD(3) * pkin(6) + t67;
t52 = qJ(5) * qJD(3) + t29;
t68 = qJD(1) * t42;
t6 = -t52 * t40 - t68;
t5 = qJD(4) * pkin(4) + t6;
t37 = t40 * qJD(1);
t7 = t52 * t42 - t37;
t47 = m(6) * (t40 * t5 - t42 * t7);
t63 = qJD(3) * t42;
t24 = -qJD(4) * mrSges(6,2) + mrSges(6,3) * t63;
t25 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t63;
t71 = t24 + t25;
t65 = qJD(3) * t40;
t22 = qJD(4) * mrSges(6,1) - mrSges(6,3) * t65;
t23 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t65;
t72 = t22 + t23;
t8 = -t29 * t40 - t68;
t9 = t29 * t42 - t37;
t87 = -t47 - m(5) * (t40 * t8 - t42 * t9) - t72 * t40 + t71 * t42;
t34 = -pkin(4) * t42 - pkin(3);
t43 = cos(qJ(3));
t66 = qJD(2) * t43;
t12 = t34 * qJD(3) + qJD(5) - t66;
t20 = (-t42 * mrSges(6,1) + t40 * mrSges(6,2)) * qJD(3);
t69 = qJD(3) * pkin(3);
t30 = -t66 - t69;
t76 = m(6) * t12;
t86 = (t20 + t76) * pkin(4) + t30 * mrSges(5,1) + t12 * mrSges(6,1) - t9 * mrSges(5,3) - t7 * mrSges(6,3) - ((Ifges(5,2) + Ifges(6,2)) * t42 + t89 * t40) * qJD(3) / 0.2e1 + (Ifges(6,6) + Ifges(5,6)) * t90;
t59 = qJD(3) * qJD(4);
t17 = (mrSges(6,1) * t40 + mrSges(6,2) * t42) * t59;
t62 = qJD(4) * t40;
t19 = (pkin(4) * t62 + t67) * qJD(3);
t81 = m(6) * t19 + t17;
t79 = -t89 * t63 / 0.2e1;
t77 = m(5) * t9;
t74 = -qJ(5) - pkin(6);
t73 = t20 + (-mrSges(5,1) * t42 + mrSges(5,2) * t40) * qJD(3);
t64 = qJD(3) * t41;
t61 = qJD(4) * t42;
t60 = qJD(5) * t42;
t58 = 0.3e1 / 0.2e1 * Ifges(6,4) + 0.3e1 / 0.2e1 * Ifges(5,4);
t57 = Ifges(5,5) / 0.2e1 + Ifges(6,5) / 0.2e1;
t56 = m(5) * pkin(6) + mrSges(5,3);
t55 = -m(5) * t8 - t23;
t54 = qJD(3) * t66;
t53 = qJD(4) * t74;
t51 = qJD(4) * t37 - t29 * t61;
t48 = (-Ifges(6,6) / 0.2e1 - Ifges(5,6) / 0.2e1) * qJD(4);
t46 = -t30 * mrSges(5,2) - t12 * mrSges(6,2) + t8 * mrSges(5,3) + t5 * mrSges(6,3) + t79 - (Ifges(5,1) + Ifges(6,1)) * t65 / 0.2e1 + (Ifges(5,5) + Ifges(6,5)) * t90;
t44 = qJD(3) ^ 2;
t28 = t42 * t54;
t27 = t74 * t42;
t26 = t74 * t40;
t18 = (mrSges(5,1) * t40 + mrSges(5,2) * t42) * t59;
t11 = -qJD(5) * t40 + t42 * t53;
t10 = t40 * t53 + t60;
t4 = -t40 * t54 + t51;
t3 = t8 * qJD(4) + t28;
t2 = (-qJ(5) * t61 + (-qJD(5) - t66) * t40) * qJD(3) + t51;
t1 = qJD(3) * t60 + t6 * qJD(4) + t28;
t13 = [m(5) * (-t3 * t40 - t4 * t42) + m(6) * (-t1 * t40 - t2 * t42) + ((mrSges(5,3) + mrSges(6,3)) * qJD(3) * (t40 ^ 2 + t42 ^ 2) - t87) * qJD(4); (-t44 * mrSges(4,2) + t87 * qJD(3) - t18 - t81) * t43 + (-t44 * mrSges(4,1) + t73 * qJD(3) + (-t71 * t40 - t72 * t42) * qJD(4) + m(5) * (qJD(3) * t30 + t3 * t42 - t4 * t40 - t8 * t61 - t9 * t62 - t54) + m(6) * (qJD(3) * t12 + t1 * t42 - t2 * t40 - t5 * t61 - t7 * t62)) * t41; -pkin(3) * t18 + t10 * t24 + t11 * t22 + t34 * t17 + m(6) * (-t1 * t27 + t7 * t10 + t5 * t11 + t19 * t34 + t2 * t26) + (t19 * mrSges(6,2) - t2 * mrSges(6,3) - t56 * t4 + (mrSges(5,2) * t64 + (m(6) * t5 + t22 - t55) * t43) * qJD(2) + (t48 + (t27 * mrSges(6,3) - t58 * t40) * qJD(3) + (-t25 - t77) * pkin(6) + t86) * qJD(4)) * t40 + (-t19 * mrSges(6,1) + t1 * mrSges(6,3) + t56 * t3 + (-mrSges(5,1) * t64 + (-m(6) * t7 - t71 - t77) * t43) * qJD(2) + (t57 * qJD(4) + t55 * pkin(6) + (-t26 * mrSges(6,3) + t58 * t42 + (0.3e1 / 0.2e1 * Ifges(5,1) + 0.3e1 / 0.2e1 * Ifges(6,1) - 0.3e1 / 0.2e1 * Ifges(6,2) - 0.3e1 / 0.2e1 * Ifges(5,2)) * t40) * qJD(3) - t46) * qJD(4)) * t42 + (-t73 - t76 + 0.2e1 * (-t30 / 0.2e1 - t69 / 0.2e1) * m(5)) * t67; t4 * mrSges(5,1) - t3 * mrSges(5,2) - t1 * mrSges(6,2) + t9 * t23 - t6 * t24 - t8 * t25 + (m(6) * pkin(4) + mrSges(6,1)) * t2 + (t22 - m(6) * (-t5 + t6)) * t7 + (((-pkin(4) * mrSges(6,3) + t57) * qJD(4) + t46 + t79) * t42 + ((Ifges(6,4) / 0.2e1 + Ifges(5,4) / 0.2e1) * t65 + t48 + (-Ifges(6,1) / 0.2e1 - Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t63 - t86) * t40) * qJD(3); (t40 * t22 - t42 * t24 + t47) * qJD(3) + t81;];
tauc = t13(:);
