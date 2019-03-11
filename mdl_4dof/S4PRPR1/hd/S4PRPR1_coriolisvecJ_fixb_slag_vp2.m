% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PRPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d4,theta1]';
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
% Datum: 2019-03-08 18:21
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PRPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:20:52
% EndTime: 2019-03-08 18:20:52
% DurationCPUTime: 0.12s
% Computational Cost: add. (117->33), mult. (228->51), div. (0->0), fcn. (74->2), ass. (0->22)
t10 = cos(qJ(4));
t8 = -qJD(2) + qJD(4);
t9 = sin(qJ(4));
t24 = (mrSges(5,1) * t9 + mrSges(5,2) * t10) * t8;
t23 = (m(4) * qJ(3) + mrSges(4,3)) * qJD(2);
t11 = -pkin(2) - pkin(3);
t7 = t11 * qJD(2) + qJD(3);
t22 = qJD(4) * t7;
t21 = t9 * qJD(3);
t20 = t10 * qJD(3);
t19 = qJ(3) * qJD(2);
t18 = qJ(3) * qJD(4);
t5 = t10 * t7 - t9 * t19;
t6 = t10 * t19 + t9 * t7;
t16 = t10 * t6 - t5 * t9;
t14 = -t9 * qJ(3) + t10 * t11;
t13 = t10 * qJ(3) + t9 * t11;
t4 = -t13 * qJD(4) - t21;
t3 = t14 * qJD(4) + t20;
t2 = -t9 * t22 + (-t10 * t18 - t21) * qJD(2);
t1 = t10 * t22 + (-t9 * t18 + t20) * qJD(2);
t12 = [0; m(5) * (t1 * t13 + t2 * t14 + t6 * t3 + t5 * t4) + 0.2e1 * qJD(3) * t23 + (-t3 * t8 + t1) * mrSges(5,2) + (t4 * t8 - t2) * mrSges(5,1); m(5) * (t16 * qJD(4) + t1 * t9 + t2 * t10) - qJD(4) * t24 + (-m(5) * t16 - t23 + t24) * qJD(2); (t5 * t8 - t1) * mrSges(5,2) + (t6 * t8 + t2) * mrSges(5,1);];
tauc  = t12(:);
