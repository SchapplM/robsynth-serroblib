% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,d3]';
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
% Datum: 2018-11-14 14:03
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRRP2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRRP2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:03:19
% EndTime: 2018-11-14 14:03:19
% DurationCPUTime: 0.21s
% Computational Cost: add. (225->45), mult. (581->69), div. (0->0), fcn. (358->4), ass. (0->28)
t37 = mrSges(4,1) + mrSges(5,1);
t25 = cos(qJ(2));
t20 = qJD(2) * pkin(2) + t25 * qJD(1);
t22 = sin(qJ(3));
t24 = cos(qJ(3));
t23 = sin(qJ(2));
t35 = qJD(1) * t23;
t16 = t22 * t20 + t24 * t35;
t30 = -t22 * t23 + t24 * t25;
t18 = t30 * qJD(1);
t33 = qJD(3) * t24;
t27 = t30 * qJD(2);
t34 = qJD(3) * t22;
t8 = t20 * t33 + (-t23 * t34 + t27) * qJD(1);
t39 = -t16 * t18 + (t16 * t33 + t22 * t8) * pkin(2);
t36 = mrSges(4,2) + mrSges(5,2);
t11 = t30 * qJD(3) + t27;
t31 = t22 * t25 + t24 * t23;
t28 = t31 * qJD(2);
t9 = -t20 * t34 + (-t23 * t33 - t28) * qJD(1);
t32 = t16 * t11 + t9 * t30 + t8 * t31;
t15 = t24 * t20 - t22 * t35;
t29 = -t36 * t8 + t37 * t9;
t21 = qJD(2) + qJD(3);
t17 = t31 * qJD(1);
t14 = t21 * pkin(3) + t15;
t12 = -t31 * qJD(3) - t28;
t1 = [(-t23 * mrSges(3,1) - t25 * mrSges(3,2)) * qJD(2) ^ 2 + m(4) * (t15 * t12 + t32) + m(5) * (t14 * t12 + t32) + (-t36 * t11 + t37 * t12) * t21; (t36 * t18 + t37 * t17 + (-t37 * t22 - t36 * t24) * qJD(3) * pkin(2)) * t21 + t29 + (t9 * (t24 * pkin(2) + pkin(3)) + (-pkin(2) * t34 + t17) * t14 + t39) * m(5) + ((-t15 * t34 + t9 * t24) * pkin(2) + t15 * t17 + t39) * m(4); (t36 * t15 + t37 * t16) * t21 + t29 + (t9 * pkin(3) + (t14 - t15) * t16) * m(5); 0;];
tauc  = t1(:);
