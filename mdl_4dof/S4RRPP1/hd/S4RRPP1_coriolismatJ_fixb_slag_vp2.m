% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RRPP1_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP1_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:51:30
% EndTime: 2018-11-14 13:51:30
% DurationCPUTime: 0.14s
% Computational Cost: add. (168->37), mult. (432->52), div. (0->0), fcn. (257->4), ass. (0->27)
t15 = sin(pkin(6));
t18 = cos(qJ(2));
t16 = cos(pkin(6));
t17 = sin(qJ(2));
t31 = t16 * t17;
t11 = (t15 * t18 + t31) * pkin(1);
t32 = t15 * t17;
t12 = (t16 * t18 - t32) * pkin(1);
t36 = -(t17 * mrSges(3,1) + t18 * mrSges(3,2)) * pkin(1) - (mrSges(4,2) - mrSges(5,3)) * t12 - (mrSges(4,1) + mrSges(5,1)) * t11;
t34 = m(5) / 0.2e1;
t33 = t15 * pkin(2);
t25 = t18 * pkin(1) + pkin(2);
t21 = pkin(1) * t31 + t15 * t25;
t20 = qJ(4) + t21;
t22 = -pkin(1) * t32 + t16 * t25;
t1 = -m(4) * (-t22 * t11 + t21 * t12) - m(5) * (t20 * t12 + (-pkin(3) - t22) * t11) - t36;
t28 = t1 * qJD(1);
t7 = m(5) * t20 + mrSges(5,3);
t27 = t7 * qJD(1);
t26 = t11 * t34;
t14 = qJ(4) + t33;
t13 = m(5) * t14 + mrSges(5,3);
t19 = -mrSges(5,3) - m(5) * (0.2e1 * qJ(4) + t11 + 0.2e1 * t33) / 0.2e1;
t4 = t26 + t19;
t24 = t4 * qJD(1) - t13 * qJD(2);
t5 = t26 - t19;
t2 = [-t1 * qJD(2) + t7 * qJD(4), t5 * qJD(4) - t28 + (m(4) * (-t16 * t11 + t15 * t12) * pkin(2) + 0.2e1 * (t14 * t12 + (-t16 * pkin(2) - pkin(3)) * t11) * t34 + t36) * qJD(2), 0, t5 * qJD(2) + t27; -t4 * qJD(4) + t28, t13 * qJD(4), 0, -t24; 0, 0, 0, 0; t4 * qJD(2) - t27, t24, 0, 0;];
Cq  = t2;
