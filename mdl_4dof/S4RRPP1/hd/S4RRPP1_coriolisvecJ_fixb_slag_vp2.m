% Calculate vector of centrifugal and coriolis load on the joints for
% S4RRPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2,theta3]';
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

function tauc = S4RRPP1_coriolisvecJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_coriolisvecJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolisvec_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:29
% EndTime: 2018-11-14 13:51:30
% DurationCPUTime: 0.22s
% Computational Cost: add. (156->52), mult. (444->73), div. (0->0), fcn. (200->4), ass. (0->37)
t38 = mrSges(4,1) + mrSges(5,1);
t36 = pkin(1) * qJD(1);
t24 = cos(pkin(6));
t26 = cos(qJ(2));
t39 = t24 * t26;
t23 = sin(pkin(6));
t25 = sin(qJ(2));
t41 = t23 * t25;
t13 = (t39 - t41) * t36;
t42 = qJD(4) - t13;
t40 = t24 * t25;
t21 = t26 * pkin(1) + pkin(2);
t37 = pkin(1) * t40 + t23 * t21;
t35 = pkin(1) * qJD(2);
t34 = pkin(1) * t41;
t33 = t25 * t36;
t32 = t23 * t33;
t19 = t35 * t39;
t31 = mrSges(3,1) * t25 + mrSges(3,2) * t26;
t22 = qJD(1) + qJD(2);
t16 = t22 * pkin(2) + t26 * t36;
t6 = t23 * t16 + t24 * t33;
t14 = -qJD(2) * t34 + t19;
t30 = t24 * t21 - t34;
t29 = pkin(1) * (t23 * t26 + t40);
t9 = qJD(1) * t19 - qJD(2) * t32;
t5 = t24 * t16 - t32;
t28 = t31 * t35;
t12 = qJD(2) * t29;
t4 = t22 * qJD(4) + t9;
t8 = qJD(1) * t12;
t27 = -t9 * mrSges(4,2) + t4 * mrSges(5,3) - qJD(1) * t28 - t38 * t8;
t11 = qJD(1) * t29;
t10 = qJD(4) + t14;
t2 = t22 * qJ(4) + t6;
t1 = -t22 * pkin(3) + qJD(4) - t5;
t3 = [m(4) * (-t5 * t12 + t6 * t14 - t8 * t30 + t9 * t37) + m(5) * (t4 * (qJ(4) + t37) + t2 * t10 + t8 * (-pkin(3) - t30) + t1 * t12) + (-t14 * mrSges(4,2) + t10 * mrSges(5,3) - t38 * t12 - t28) * t22 + t27; (t13 * mrSges(4,2) + t42 * mrSges(5,3) + t38 * t11 + t31 * t36) * t22 + t27 + (t4 * (t23 * pkin(2) + qJ(4)) + t8 * (-t24 * pkin(2) - pkin(3)) - t1 * t11 + t42 * t2) * m(5) + ((t23 * t9 - t8 * t24) * pkin(2) + t5 * t11 - t6 * t13) * m(4); 0; -t22 ^ 2 * mrSges(5,3) + (-t2 * t22 + t8) * m(5);];
tauc  = t3(:);
