% Calculate vector of centrifugal and coriolis load on the joints for
% S4PRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d2,theta3]';
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
% Datum: 2018-11-14 14:09
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PRPP4_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PRPP4_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:09:18
% EndTime: 2018-11-14 14:09:19
% DurationCPUTime: 0.15s
% Computational Cost: add. (98->42), mult. (295->60), div. (0->0), fcn. (178->4), ass. (0->27)
t17 = sin(pkin(5));
t18 = cos(pkin(5));
t19 = sin(qJ(2));
t20 = cos(qJ(2));
t12 = t17 * t19 - t18 * t20;
t10 = t12 * qJD(1);
t27 = qJD(4) + t10;
t13 = t17 * t20 + t18 * t19;
t9 = t13 * qJD(2);
t6 = qJD(1) * t9;
t26 = t6 * t12;
t25 = mrSges(4,1) + mrSges(5,1);
t24 = qJD(1) * t19;
t23 = t20 * qJD(1);
t22 = t17 * t24;
t16 = qJD(2) * pkin(2) + t23;
t5 = t17 * t16 + t18 * t24;
t4 = t18 * t16 - t22;
t21 = qJD(2) ^ 2;
t15 = t18 * qJD(2) * t23;
t11 = t12 * qJD(2);
t8 = t13 * qJD(1);
t7 = -qJD(2) * t22 + t15;
t3 = t15 + (qJD(4) - t22) * qJD(2);
t2 = qJD(2) * qJ(4) + t5;
t1 = -qJD(2) * pkin(3) + qJD(4) - t4;
t14 = [(-t19 * mrSges(3,1) - t20 * mrSges(3,2)) * t21 + m(4) * (-t5 * t11 + t7 * t13 - t4 * t9 + t26) + m(5) * (t1 * t9 - t2 * t11 + t3 * t13 + t26) + (-t25 * t9 - (-mrSges(4,2) + mrSges(5,3)) * t11) * qJD(2); -t7 * mrSges(4,2) + t3 * mrSges(5,3) - t25 * t6 + (-t10 * mrSges(4,2) + t27 * mrSges(5,3) + t25 * t8) * qJD(2) + (t3 * (t17 * pkin(2) + qJ(4)) + t6 * (-t18 * pkin(2) - pkin(3)) - t1 * t8 + t27 * t2) * m(5) + ((t17 * t7 - t6 * t18) * pkin(2) + t5 * t10 + t4 * t8) * m(4); 0; -t21 * mrSges(5,3) + (-t2 * qJD(2) + t6) * m(5);];
tauc  = t14(:);
