% Calculate vector of centrifugal and Coriolis load on the joints for
% S4PPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta1]';
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
% Datum: 2019-03-08 18:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function tauc = S4PPRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:12:14
% EndTime: 2019-03-08 18:12:14
% DurationCPUTime: 0.08s
% Computational Cost: add. (27->19), mult. (78->25), div. (0->0), fcn. (23->2), ass. (0->12)
t4 = sin(qJ(3));
t10 = qJD(2) * t4;
t3 = qJD(3) * qJ(4) + t10;
t5 = cos(qJ(3));
t11 = t3 * t5;
t9 = qJD(2) * t5;
t8 = qJD(3) * pkin(3);
t7 = qJD(4) - t9;
t6 = qJD(3) ^ 2;
t2 = t7 - t8;
t1 = (qJD(4) + t9) * qJD(3);
t12 = [0; m(5) * (t1 * t4 + (t11 + (t2 - t9) * t4) * qJD(3)) + ((-mrSges(4,2) + mrSges(5,3)) * t5 + (-mrSges(4,1) - mrSges(5,1)) * t4) * t6; (t7 * qJD(3) + t1) * mrSges(5,3) + (qJ(4) * t1 + qJD(4) * t3 - t8 * t10 - (t2 * t4 + t11) * qJD(2)) * m(5); -t6 * mrSges(5,3) + (-t3 + t10) * qJD(3) * m(5);];
tauc  = t12(:);
