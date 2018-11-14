% Calculate vector of centrifugal and coriolis load on the joints for
% S4PPRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d3,theta2]';
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
% Datum: 2018-11-14 14:08
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function tauc = S4PPRP5_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRP5_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 14:07:13
% EndTime: 2018-11-14 14:07:13
% DurationCPUTime: 0.12s
% Computational Cost: add. (80->30), mult. (257->43), div. (0->0), fcn. (164->4), ass. (0->23)
t14 = cos(pkin(5));
t16 = cos(qJ(3));
t19 = t16 * t14;
t13 = sin(pkin(5));
t15 = sin(qJ(3));
t20 = t15 * t13;
t10 = -t19 + t20;
t6 = t10 * qJD(1);
t22 = qJD(4) + t6;
t11 = t16 * t13 + t15 * t14;
t9 = t11 * qJD(3);
t5 = qJD(1) * t9;
t21 = t5 * t10;
t18 = mrSges(4,1) + mrSges(5,1);
t17 = qJD(1) * t20;
t7 = t11 * qJD(1);
t12 = qJD(3) * qJD(1) * t19;
t8 = t10 * qJD(3);
t4 = -qJD(3) * t17 + t12;
t3 = qJD(3) * qJ(4) + t7;
t2 = -qJD(3) * pkin(3) + t22;
t1 = t12 + (qJD(4) - t17) * qJD(3);
t23 = [m(4) * (t4 * t11 + t6 * t9 - t7 * t8 + t21) + m(5) * (t1 * t11 + t2 * t9 - t3 * t8 + t21) + (-t18 * t9 - (-mrSges(4,2) + mrSges(5,3)) * t8) * qJD(3); 0; -t4 * mrSges(4,2) + t1 * mrSges(5,3) - t18 * t5 + (-t6 * mrSges(4,2) + t22 * mrSges(5,3) + t18 * t7) * qJD(3) + (-t5 * pkin(3) + t1 * qJ(4) - t2 * t7 + t22 * t3) * m(5); -qJD(3) ^ 2 * mrSges(5,3) + (-t3 * qJD(3) + t5) * m(5);];
tauc  = t23(:);
