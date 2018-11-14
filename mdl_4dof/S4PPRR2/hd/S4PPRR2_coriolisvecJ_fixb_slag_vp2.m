% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,d4,theta2]';
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
% Datum: 2018-11-14 14:00
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRR2_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:59:23
% EndTime: 2018-11-14 13:59:23
% DurationCPUTime: 0.17s
% Computational Cost: add. (217->47), mult. (636->75), div. (0->0), fcn. (498->6), ass. (0->31)
t21 = sin(pkin(6));
t22 = cos(pkin(6));
t24 = sin(qJ(3));
t26 = cos(qJ(3));
t18 = -t24 * t21 + t26 * t22;
t20 = qJD(3) + qJD(4);
t33 = t20 * mrSges(5,1);
t17 = t26 * t21 + t24 * t22;
t14 = t17 * qJD(1);
t23 = sin(qJ(4));
t32 = t23 * t14;
t25 = cos(qJ(4));
t30 = t25 * t14;
t13 = t18 * qJD(1);
t10 = qJD(3) * pkin(3) + t13;
t6 = t25 * t10 - t32;
t7 = t23 * t10 + t30;
t28 = t25 * t17 + t23 * t18;
t27 = -t23 * t17 + t25 * t18;
t16 = t17 * qJD(3);
t15 = t18 * qJD(3);
t12 = qJD(1) * t16;
t11 = qJD(1) * t15;
t9 = t25 * t13 - t32;
t8 = -t23 * t13 - t30;
t5 = -t28 * qJD(4) - t23 * t15 - t25 * t16;
t4 = t27 * qJD(4) + t25 * t15 - t23 * t16;
t3 = -qJD(4) * t7 - t23 * t11 - t25 * t12;
t2 = qJD(4) * t6 + t25 * t11 - t23 * t12;
t1 = t3 * mrSges(5,1);
t19 = [(t5 * mrSges(5,1) - t4 * mrSges(5,2)) * t20 + (-t16 * mrSges(4,1) - t15 * mrSges(4,2)) * qJD(3) + m(4) * (t11 * t17 - t12 * t18 - t13 * t16 + t14 * t15) + m(5) * (t2 * t28 + t3 * t27 + t7 * t4 + t6 * t5); 0; t1 - t8 * t33 - m(5) * (t6 * t8 + t7 * t9) + (t9 * t20 - t2) * mrSges(5,2) + (t13 * qJD(3) - t11) * mrSges(4,2) + (t14 * qJD(3) - t12) * mrSges(4,1) + (m(5) * (t2 * t23 + t3 * t25 + (-t23 * t6 + t25 * t7) * qJD(4)) + (-mrSges(5,1) * t23 - mrSges(5,2) * t25) * qJD(4) * t20) * pkin(3); t7 * t33 + t1 + (t20 * t6 - t2) * mrSges(5,2);];
tauc  = t19(:);
