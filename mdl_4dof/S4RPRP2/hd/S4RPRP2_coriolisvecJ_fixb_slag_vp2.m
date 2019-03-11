% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3]';
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
% Datum: 2019-03-08 18:31
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:30:42
% EndTime: 2019-03-08 18:30:42
% DurationCPUTime: 0.20s
% Computational Cost: add. (247->45), mult. (459->63), div. (0->0), fcn. (151->2), ass. (0->28)
t42 = qJD(1) - qJD(3);
t35 = mrSges(4,2) + mrSges(5,2);
t21 = -pkin(1) - pkin(2);
t16 = t21 * qJD(1) + qJD(2);
t19 = sin(qJ(3));
t27 = qJ(2) * qJD(3);
t20 = cos(qJ(3));
t29 = t20 * qJD(2);
t31 = qJD(3) * t20;
t6 = t16 * t31 + (-t19 * t27 + t29) * qJD(1);
t41 = t35 * t6;
t28 = qJ(2) * qJD(1);
t14 = t19 * t16 + t20 * t28;
t30 = t19 * qJD(2);
t32 = qJD(3) * t19;
t7 = -t16 * t32 + (-t20 * t27 - t30) * qJD(1);
t40 = t6 * t19 + t7 * t20 + (-qJD(1) * t20 + t31) * t14;
t24 = -t19 * qJ(2) + t20 * t21;
t39 = qJD(1) * t19 - t32;
t11 = t24 * qJD(3) + t29;
t23 = t20 * qJ(2) + t19 * t21;
t37 = t14 * t11 + t6 * t23;
t36 = mrSges(4,1) + mrSges(5,1);
t25 = m(3) * qJ(2) + mrSges(3,3);
t13 = t20 * t16 - t19 * t28;
t12 = -t23 * qJD(3) - t30;
t8 = -pkin(3) * t42 + t13;
t1 = [-t36 * t7 + t41 + 0.2e1 * t25 * qJD(2) * qJD(1) + m(4) * (t13 * t12 + t7 * t24 + t37) + m(5) * (t7 * (-pkin(3) + t24) + t8 * t12 + t37) - (-t35 * t11 + t36 * t12) * t42; -t25 * qJD(1) ^ 2 + (t39 * t8 + t40) * m(5) + (t39 * t13 + t40) * m(4) - t42 ^ 2 * (t36 * t19 + t35 * t20); -t41 + (m(5) * pkin(3) + t36) * t7 - t35 * t13 * t42 + (-m(5) * (t13 - t8) - t36 * t42) * t14; 0;];
tauc  = t1(:);
