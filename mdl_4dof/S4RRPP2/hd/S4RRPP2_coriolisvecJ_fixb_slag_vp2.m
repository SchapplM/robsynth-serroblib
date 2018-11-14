% Calculate vector of centrifugal and coriolis load on the joints for
% S4RRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RRPP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:24
% EndTime: 2018-11-14 13:52:25
% DurationCPUTime: 0.19s
% Computational Cost: add. (124->41), mult. (254->48), div. (0->0), fcn. (58->2), ass. (0->28)
t15 = cos(qJ(2));
t28 = pkin(1) * qJD(1);
t16 = -t15 * t28 + qJD(3);
t29 = mrSges(5,2) + mrSges(4,3);
t13 = qJD(1) + qJD(2);
t14 = sin(qJ(2));
t23 = t14 * t28;
t11 = t13 * qJ(3) + t23;
t27 = pkin(1) * qJD(2);
t20 = qJD(1) * t27;
t17 = t15 * t20;
t9 = t13 * qJD(3) + t17;
t38 = t9 * qJ(3) + t16 * t11;
t26 = mrSges(3,1) + mrSges(4,1) + mrSges(5,1);
t37 = t26 * t14;
t36 = t29 * t13;
t35 = m(4) + m(5);
t34 = -pkin(2) - pkin(3);
t24 = t15 * t27;
t12 = qJD(3) + t24;
t33 = t9 * (t14 * pkin(1) + qJ(3)) + t11 * t12;
t31 = t29 * t9;
t25 = t14 * t27;
t21 = -t15 * pkin(1) - pkin(2);
t18 = t14 * t20;
t10 = -t13 * pkin(2) + t16;
t3 = t34 * t13 + t16;
t1 = [m(4) * ((t21 * qJD(1) + t10) * t25 + t33) + m(5) * ((t3 + (-pkin(3) + t21) * qJD(1)) * t25 + t33) + t31 + t12 * t36 + (-t13 * t24 - t17) * mrSges(3,2) + (-t13 * t25 - t18) * t26; (-mrSges(3,2) * t15 - t37) * t20 + (t29 * qJD(3) + ((mrSges(3,2) - t29) * t15 + t37) * t28) * t13 + t31 + (t34 * t18 - t3 * t23 + t38) * m(5) + (-pkin(2) * t18 - t10 * t23 + t38) * m(4); t35 * t18 + (-t35 * t11 - t36) * t13; 0;];
tauc  = t1(:);
