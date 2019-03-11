% Calculate vector of centrifugal and Coriolis load on the joints for
% S4RPPR2
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
% Datum: 2019-03-08 18:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4RPPR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:28:26
% EndTime: 2019-03-08 18:28:27
% DurationCPUTime: 0.24s
% Computational Cost: add. (289->45), mult. (546->68), div. (0->0), fcn. (266->4), ass. (0->32)
t44 = qJD(1) - qJD(4);
t20 = sin(pkin(6));
t21 = cos(pkin(6));
t43 = qJD(1) * (m(3) * qJ(2) + t20 * mrSges(4,1) + t21 * mrSges(4,2) + mrSges(3,3));
t22 = sin(qJ(4));
t23 = cos(qJ(4));
t29 = t23 * t20 + t22 * t21;
t42 = t44 * t29;
t28 = -t22 * t20 + t23 * t21;
t41 = t44 * t28;
t24 = -pkin(1) - pkin(2);
t36 = t20 * qJ(2);
t35 = qJD(1) * qJ(2);
t33 = t21 * t24 - t36;
t17 = t24 * qJD(1) + qJD(2);
t16 = t21 * t17;
t9 = t20 * t17 + t21 * t35;
t32 = -(-t20 * t35 + t16) * t20 + t9 * t21;
t7 = t16 + (-pkin(3) - t36) * qJD(1);
t5 = -t22 * t9 + t23 * t7;
t6 = t22 * t7 + t23 * t9;
t14 = -pkin(3) + t33;
t15 = t21 * qJ(2) + t20 * t24;
t31 = t23 * t14 - t22 * t15;
t30 = t22 * t14 + t23 * t15;
t27 = t29 * qJD(2);
t26 = t28 * qJD(2);
t4 = -t30 * qJD(4) - t27;
t3 = t31 * qJD(4) + t26;
t2 = -qJD(1) * t27 - t6 * qJD(4);
t1 = qJD(1) * t26 + t5 * qJD(4);
t8 = [m(5) * (t1 * t30 + t2 * t31 + t6 * t3 + t5 * t4) + (t3 * t44 + t1) * mrSges(5,2) + (-t4 * t44 - t2) * mrSges(5,1) + (m(4) * ((t21 * t15 - t20 * t33) * qJD(1) + t32) + 0.2e1 * t43) * qJD(2); -(t42 * mrSges(5,1) + t41 * mrSges(5,2)) * t44 + (-m(4) * t32 - t43) * qJD(1) + (t1 * t29 + t2 * t28 - t41 * t6 + t42 * t5) * m(5); 0; (-t44 * t5 - t1) * mrSges(5,2) + (-t44 * t6 + t2) * mrSges(5,1);];
tauc  = t8(:);
