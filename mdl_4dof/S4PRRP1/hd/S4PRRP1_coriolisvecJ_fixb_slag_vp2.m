% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3,theta1]';
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
% Datum: 2018-11-14 13:43
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRRP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:43:27
% EndTime: 2018-11-14 13:43:27
% DurationCPUTime: 0.13s
% Computational Cost: add. (63->31), mult. (142->45), div. (0->0), fcn. (33->2), ass. (0->18)
t18 = mrSges(4,1) + mrSges(5,1);
t17 = pkin(2) * qJD(2);
t16 = pkin(2) * qJD(3);
t6 = qJD(2) + qJD(3);
t15 = t6 * qJD(4);
t7 = sin(qJ(3));
t14 = t7 * t16;
t8 = cos(qJ(3));
t13 = t8 * t16;
t12 = qJD(2) * t16;
t10 = t7 * t12;
t9 = t8 * t12;
t5 = qJD(4) + t13;
t4 = qJ(4) * t6 + t7 * t17;
t3 = -t6 * pkin(3) - t8 * t17 + qJD(4);
t2 = t9 + t15;
t1 = t2 * mrSges(5,3);
t11 = [0; t1 + t5 * t6 * mrSges(5,3) + m(5) * (t2 * (pkin(2) * t7 + qJ(4)) + t4 * t5 + (t3 + (-t8 * pkin(2) - pkin(3)) * qJD(2)) * t14) + (-t6 * t13 - t9) * mrSges(4,2) + t18 * (-t6 * t14 - t10); mrSges(5,3) * t15 + t1 + m(5) * (t2 * qJ(4) + qJD(4) * t4) + (-m(5) * (t3 * t7 + t4 * t8) + ((mrSges(4,2) - mrSges(5,3)) * t8 + t18 * t7) * t6 + (-t8 * mrSges(4,2) + (-m(5) * pkin(3) - t18) * t7) * qJD(3)) * t17; -t6 ^ 2 * mrSges(5,3) + (-t4 * t6 + t10) * m(5);];
tauc  = t11(:);
