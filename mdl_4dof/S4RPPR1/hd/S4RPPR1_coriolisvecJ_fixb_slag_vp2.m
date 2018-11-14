% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d4,theta2]';
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
% Datum: 2018-11-14 13:47
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPPR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:46:34
% EndTime: 2018-11-14 13:46:34
% DurationCPUTime: 0.13s
% Computational Cost: add. (173->33), mult. (329->50), div. (0->0), fcn. (129->4), ass. (0->20)
t11 = -qJD(1) + qJD(4);
t14 = sin(qJ(4));
t15 = cos(qJ(4));
t24 = (mrSges(5,1) * t14 + mrSges(5,2) * t15) * t11;
t10 = sin(pkin(6)) * pkin(1) + qJ(3);
t8 = qJD(1) * t10;
t23 = m(4) * t8 + qJD(1) * mrSges(4,3);
t20 = qJD(1) * qJD(3);
t9 = -cos(pkin(6)) * pkin(1) - pkin(2) - pkin(3);
t7 = t9 * qJD(1) + qJD(3);
t3 = -t14 * t8 + t15 * t7;
t4 = t14 * t7 + t15 * t8;
t19 = -t14 * t3 + t15 * t4;
t18 = t15 * t10 + t14 * t9;
t17 = -t14 * t10 + t15 * t9;
t6 = -t14 * qJD(3) - t18 * qJD(4);
t5 = t15 * qJD(3) + t17 * qJD(4);
t2 = -qJD(4) * t4 - t14 * t20;
t1 = qJD(4) * t3 + t15 * t20;
t12 = [m(5) * (t1 * t18 + t2 * t17 + t3 * t6 + t4 * t5) + (-t5 * t11 + t1) * mrSges(5,2) + (t6 * t11 - t2) * mrSges(5,1) + 0.2e1 * t23 * qJD(3); 0; m(5) * (t19 * qJD(4) + t1 * t14 + t2 * t15) - qJD(4) * t24 + (-m(5) * t19 - t23 + t24) * qJD(1); (t11 * t3 - t1) * mrSges(5,2) + (t11 * t4 + t2) * mrSges(5,1);];
tauc  = t12(:);
