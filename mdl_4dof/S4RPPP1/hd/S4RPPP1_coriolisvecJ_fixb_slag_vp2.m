% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:25
% EndTime: 2018-11-14 13:45:26
% DurationCPUTime: 0.50s
% Computational Cost: add. (183->84), mult. (649->132), div. (0->0), fcn. (393->4), ass. (0->52)
t32 = cos(pkin(4));
t30 = sin(pkin(4));
t31 = cos(pkin(6));
t59 = t30 * t31;
t65 = (-mrSges(4,1) - mrSges(5,1)) * t59 + (-mrSges(5,2) - mrSges(4,3)) * t32;
t50 = t30 * qJD(1);
t29 = sin(pkin(6));
t9 = (mrSges(4,2) * t31 - mrSges(4,3) * t29) * t50;
t45 = -pkin(1) * t31 - pkin(2);
t64 = (-qJ(4) + t45) * t32;
t58 = t65 * qJD(1);
t60 = t29 * t30;
t63 = ((mrSges(4,1) + mrSges(3,3)) * t60 + (-mrSges(3,1) + mrSges(4,2)) * t32) * qJD(1);
t55 = t32 * qJ(3);
t61 = pkin(1) * t32;
t44 = qJ(2) * t50;
t54 = qJD(1) * t32;
t48 = pkin(1) * t54;
t11 = t29 * t48 + t31 * t44;
t56 = qJ(2) * t30;
t57 = t29 * t61 + t31 * t56;
t53 = qJD(2) * t30;
t52 = t29 * qJD(2);
t51 = t29 * qJD(3);
t23 = t29 * t44;
t49 = qJD(3) + t23;
t47 = t30 * t52;
t46 = t31 * t53;
t43 = -qJ(3) * t29 - pkin(1);
t41 = t45 * t32;
t40 = t11 * t31 - (t31 * t48 - t23) * t29;
t37 = mrSges(5,1) * t60 - t32 * mrSges(5,3);
t36 = (-mrSges(5,2) * t29 - mrSges(5,3) * t31) * t30;
t22 = t32 * qJD(3) + t46;
t21 = -t32 * qJD(4) + t47;
t12 = (-qJD(4) * t31 - t51) * t30;
t34 = (-pkin(2) * t31 + t43) * t30;
t33 = (-pkin(2) - qJ(4)) * t31 + t43;
t16 = (-t32 * mrSges(3,2) + mrSges(3,3) * t59) * qJD(1);
t26 = t29 * t56;
t17 = t37 * qJD(1);
t14 = t22 * qJD(1);
t13 = t21 * qJD(1);
t10 = qJD(1) * t36;
t7 = qJD(1) * t12;
t6 = -qJ(3) * t54 - t11;
t5 = qJD(1) * t34 + qJD(2);
t4 = qJD(1) * t41 + t49;
t3 = t33 * t50 + qJD(2);
t2 = qJD(4) + (pkin(3) * t59 + t55) * qJD(1) + t11;
t1 = (pkin(3) * t60 + t64) * qJD(1) + t49;
t8 = [t12 * t10 + t21 * t17 + 0.2e1 * t16 * t46 + t7 * t36 + m(3) * ((t31 * t57 - t29 * (t31 * t61 - t26)) * qJD(1) + t40) * t53 + m(5) * (t1 * t21 + t3 * t12) + 0.2e1 * t63 * t47 + (t37 + m(5) * (t26 + t64)) * t13 + ((-m(4) - m(5)) * (-t55 - t57) - t65) * t14 + (-0.2e1 * t9 * t51 + m(4) * ((t4 * qJD(2) - t5 * qJD(3)) * t29 + ((t26 + t41) * t52 - t34 * t51) * qJD(1)) + m(5) * (t7 * t33 + (t13 * t29 + t14 * t31) * pkin(3))) * t30 + (-m(4) * t6 + m(5) * t2 - t58) * t22; m(5) * t7 + ((-t16 + t58) * t31 + (-m(4) * qJD(3) - t17 - t63) * t29 - m(3) * t40 - m(4) * (t4 * t29 - t31 * t6) - m(5) * (t1 * t29 + t2 * t31)) * t50; m(5) * t13 + (t58 * t32 + (m(4) * qJD(2) + t10 + t9) * t60 - m(4) * (-t32 * t6 - t5 * t60) - m(5) * (t2 * t32 - t3 * t60)) * qJD(1); m(5) * t14 + (t32 * t17 + t10 * t59 - m(5) * (-t1 * t32 - t3 * t59)) * qJD(1);];
tauc  = t8(:);
