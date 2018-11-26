% Calculate vector of centrifugal and coriolis load on the joints for
% S2RR1
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
%   joint torques required to compensate coriolis and centrifugal load

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-16 16:44
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S2RR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(2,1),zeros(2,1),zeros(1,1),zeros(3,1),zeros(3,3),zeros(3,6)}
assert(isreal(qJ) && all(size(qJ) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [2x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [2 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [2x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [1x1] (double)');
assert( isreal(m) && all(size(m) == [3 1]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: m has to be [3x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [3,3]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [3x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [3 6]), ...
  'S2RR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [3x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-16 16:44:19
% EndTime: 2018-11-16 16:44:20
% DurationCPUTime: 0.10s
% Computational Cost: add. (42->26), mult. (138->50), div. (0->0), fcn. (52->2), ass. (0->14)
t5 = sin(qJ(2));
t13 = t5 / 0.2e1;
t6 = cos(qJ(2));
t12 = t6 / 0.2e1;
t11 = Ifges(3,4) * t5;
t10 = Ifges(3,4) * t6;
t9 = mrSges(3,3) * qJD(1);
t8 = Ifges(3,5) * qJD(2);
t7 = Ifges(3,6) * qJD(2);
t4 = -qJD(2) * mrSges(3,2) - t6 * t9;
t3 = qJD(2) * mrSges(3,1) + t5 * t9;
t2 = t8 + (-Ifges(3,1) * t5 - t10) * qJD(1);
t1 = t7 + (-Ifges(3,2) * t6 - t11) * qJD(1);
t14 = [((-t2 / 0.2e1 + pkin(1) * t3 - t8 / 0.2e1 + 0.3e1 / 0.2e1 * qJD(1) * t10) * t6 + (t1 / 0.2e1 + pkin(1) * t4 + t7 / 0.2e1 + (-0.3e1 / 0.2e1 * t11 + (0.3e1 / 0.2e1 * Ifges(3,1) - 0.3e1 / 0.2e1 * Ifges(3,2)) * t6) * qJD(1)) * t5) * qJD(2); (t2 * t12 - t5 * t1 / 0.2e1 + ((-Ifges(3,1) * t6 + t11) * t13 + (Ifges(3,2) * t5 - t10) * t12) * qJD(1) + (-t6 * t3 - t5 * t4) * pkin(1) + (-Ifges(3,5) * t6 / 0.2e1 + Ifges(3,6) * t13 + (t6 * mrSges(3,1) - t5 * mrSges(3,2)) * pkin(1)) * qJD(2)) * qJD(1);];
tauc  = t14(:);
