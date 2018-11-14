% Calculate vector of centrifugal and coriolis load on the joints for
% S4RPRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2018-11-14 13:51
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4RPRR1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:50:29
% EndTime: 2018-11-14 13:50:29
% DurationCPUTime: 0.21s
% Computational Cost: add. (391->55), mult. (1016->82), div. (0->0), fcn. (546->6), ass. (0->37)
t23 = cos(pkin(7)) * pkin(1) + pkin(2);
t21 = t23 * qJD(1);
t29 = sin(qJ(3));
t31 = cos(qJ(3));
t43 = pkin(1) * sin(pkin(7));
t36 = qJD(1) * t43;
t14 = t31 * t21 - t29 * t36;
t34 = t31 * t23 - t29 * t43;
t12 = t14 * qJD(3);
t15 = t29 * t21 + t31 * t36;
t13 = t15 * qJD(3);
t28 = sin(qJ(4));
t30 = cos(qJ(4));
t25 = qJD(1) + qJD(3);
t10 = t25 * pkin(3) + t14;
t38 = t30 * t15;
t7 = t28 * t10 + t38;
t3 = -qJD(4) * t7 - t28 * t12 - t30 * t13;
t1 = t3 * mrSges(5,1);
t42 = -t13 * mrSges(4,1) + t1;
t24 = qJD(4) + t25;
t41 = t24 * mrSges(5,1);
t40 = t25 * mrSges(4,1);
t39 = t28 * t15;
t6 = t30 * t10 - t39;
t18 = pkin(3) + t34;
t19 = t29 * t23 + t31 * t43;
t33 = t30 * t18 - t28 * t19;
t32 = t28 * t18 + t30 * t19;
t17 = t19 * qJD(3);
t16 = t34 * qJD(3);
t9 = t30 * t14 - t39;
t8 = -t28 * t14 - t38;
t5 = -t32 * qJD(4) - t28 * t16 - t30 * t17;
t4 = t33 * qJD(4) + t30 * t16 - t28 * t17;
t2 = qJD(4) * t6 + t30 * t12 - t28 * t13;
t11 = [-t17 * t40 + t5 * t41 + (-t4 * t24 - t2) * mrSges(5,2) + (-t16 * t25 - t12) * mrSges(4,2) + m(4) * (t12 * t19 - t13 * t34 - t14 * t17 + t15 * t16) + m(5) * (t2 * t32 + t3 * t33 + t7 * t4 + t6 * t5) + t42; 0; t15 * t40 - t8 * t41 - m(5) * (t6 * t8 + t7 * t9) + (t9 * t24 - t2) * mrSges(5,2) + (t14 * t25 - t12) * mrSges(4,2) + (m(5) * (t2 * t28 + t3 * t30 + (-t28 * t6 + t30 * t7) * qJD(4)) + (-mrSges(5,1) * t28 - mrSges(5,2) * t30) * qJD(4) * t24) * pkin(3) + t42; t7 * t41 + t1 + (t24 * t6 - t2) * mrSges(5,2);];
tauc  = t11(:);
