% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:48:07
% EndTime: 2019-12-31 16:48:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (534->97), mult. (1225->152), div. (0->0), fcn. (596->6), ass. (0->52)
t42 = sin(qJ(4));
t44 = cos(qJ(4));
t75 = t44 * mrSges(5,1) - t42 * mrSges(5,2) + mrSges(4,1);
t38 = cos(pkin(7)) * pkin(1) + pkin(2);
t43 = sin(qJ(3));
t45 = cos(qJ(3));
t74 = pkin(1) * sin(pkin(7));
t67 = t43 * t38 + t45 * t74;
t34 = t38 * qJD(1);
t59 = qJD(1) * t74;
t17 = t45 * t34 - t43 * t59;
t13 = t17 * qJD(3);
t18 = t43 * t34 + t45 * t59;
t39 = qJD(1) + qJD(3);
t9 = t39 * pkin(6) + t18;
t5 = t44 * qJD(2) - t42 * t9;
t66 = qJD(4) * t5;
t2 = t44 * t13 + t66;
t73 = t2 * t44;
t6 = t42 * qJD(2) + t44 * t9;
t65 = qJD(4) * t6;
t3 = -t42 * t13 - t65;
t72 = t3 * t42;
t71 = Ifges(5,4) * t42;
t69 = t39 * t42;
t68 = t39 * t44;
t64 = Ifges(5,5) * qJD(4);
t63 = Ifges(5,6) * qJD(4);
t29 = pkin(6) + t67;
t62 = qJD(4) * t29;
t61 = qJD(4) * t42;
t60 = qJD(4) * t44;
t57 = t60 / 0.2e1;
t56 = t75 * t39;
t54 = -t42 * t5 + t44 * t6;
t32 = qJD(4) * mrSges(5,1) - mrSges(5,3) * t69;
t33 = -qJD(4) * mrSges(5,2) + mrSges(5,3) * t68;
t52 = -t42 * t32 + t44 * t33;
t51 = t45 * t38 - t43 * t74;
t50 = (Ifges(5,2) * t44 + t71) * t39;
t49 = (mrSges(5,1) * t42 + mrSges(5,2) * t44) * qJD(4);
t14 = t18 * qJD(3);
t23 = t50 + t63;
t36 = Ifges(5,4) * t68;
t24 = Ifges(5,1) * t69 + t36 + t64;
t8 = -t39 * pkin(3) - t17;
t46 = mrSges(5,3) * t73 + t24 * t57 + t8 * t49 + qJD(4) ^ 2 * (Ifges(5,5) * t44 - Ifges(5,6) * t42) / 0.2e1 - (t23 + t50) * t61 / 0.2e1 - t75 * t14 + ((Ifges(5,1) * t44 - t71) * t61 + (0.3e1 * Ifges(5,4) * t44 + (Ifges(5,1) - 0.2e1 * Ifges(5,2)) * t42) * t57) * t39;
t28 = -pkin(3) - t51;
t25 = t39 * t49;
t22 = t67 * qJD(3);
t21 = t51 * qJD(3);
t1 = [t46 + (-mrSges(5,3) * t66 - t32 * t62 + t21 * t33 + m(5) * (t2 * t29 + t21 * t6 - t5 * t62)) * t44 + (-t33 * t62 - t21 * t32 + m(5) * (-t21 * t5 - t29 * t3 - t6 * t62) + (-t3 - t65) * mrSges(5,3)) * t42 - t56 * t22 + m(5) * (t14 * t28 + t8 * t22) + m(4) * (t13 * t67 - t14 * t51 - t17 * t22 + t18 * t21) + (-t21 * t39 - t13) * mrSges(4,2) + t28 * t25; m(5) * (t2 * t42 + t3 * t44) + (m(5) * t54 + (-t42 ^ 2 - t44 ^ 2) * t39 * mrSges(5,3) + t52) * qJD(4); t46 + (-t32 * t60 - t33 * t61) * pkin(6) + (-t72 + (-t6 * t42 - t5 * t44) * qJD(4)) * mrSges(5,3) + t56 * t18 + (t39 * mrSges(4,2) - t52) * t17 - t13 * mrSges(4,2) - pkin(3) * t25 + ((-t5 * t60 - t6 * t61 - t72 + t73) * pkin(6) - t54 * t17 - t8 * t18 - t14 * pkin(3)) * m(5); t3 * mrSges(5,1) - t2 * mrSges(5,2) + t6 * t32 - t5 * t33 + ((t64 / 0.2e1 - t8 * mrSges(5,2) - t24 / 0.2e1 - t36 / 0.2e1 + t5 * mrSges(5,3)) * t44 + (-t63 / 0.2e1 - t8 * mrSges(5,1) + t23 / 0.2e1 + t6 * mrSges(5,3) + (t71 / 0.2e1 + (-Ifges(5,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t44) * t39) * t42) * t39;];
tauc = t1(:);
