% Calculate vector of centrifugal and Coriolis load on the joints for
% S2RR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [2x1]
%   Generalized joint coordinates (joint angles)
% qJD [2x1]
%   Generalized joint velocities
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[d2]';
% m_mdh [3x1]
%   mass of all robot links (including the base)
% mrSges [3x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [3x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% tauc [2x1]
%   joint torques required to compensate Coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:01
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S2RR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [3 1]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:00:52
% EndTime: 2019-03-08 18:00:53
% DurationCPUTime: 0.10s
% Computational Cost: add. (42->26), mult. (138->48), div. (0->0), fcn. (52->2), ass. (0->15)
t6 = sin(qJ(2));
t15 = -t6 / 0.2e1;
t14 = t6 / 0.2e1;
t12 = Ifges(3,4) * t6;
t7 = cos(qJ(2));
t11 = Ifges(3,2) * t7;
t10 = mrSges(3,3) * qJD(1);
t9 = Ifges(3,5) * qJD(2);
t8 = Ifges(3,6) * qJD(2);
t5 = Ifges(3,4) * t7 * qJD(1);
t4 = -qJD(2) * mrSges(3,2) + t7 * t10;
t3 = qJD(2) * mrSges(3,1) - t6 * t10;
t2 = Ifges(3,1) * qJD(1) * t6 + t5 + t9;
t1 = t8 + (t11 + t12) * qJD(1);
t13 = [((t2 / 0.2e1 - pkin(1) * t3 + t9 / 0.2e1 + 0.3e1 / 0.2e1 * t5) * t7 + (-t1 / 0.2e1 - pkin(1) * t4 - t8 / 0.2e1 + (-0.3e1 / 0.2e1 * t12 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t7) * qJD(1)) * t6) * qJD(2); (t1 * t14 + ((Ifges(3,1) * t7 - t12) * t15 + t11 * t14) * qJD(1) + (t7 * t3 + t6 * t4) * pkin(1) + (Ifges(3,5) * t7 / 0.2e1 + Ifges(3,6) * t15 + (-t7 * mrSges(3,1) + t6 * mrSges(3,2)) * pkin(1)) * qJD(2) - (t2 + t5) * t7 / 0.2e1) * qJD(1);];
tauc  = t13(:);
