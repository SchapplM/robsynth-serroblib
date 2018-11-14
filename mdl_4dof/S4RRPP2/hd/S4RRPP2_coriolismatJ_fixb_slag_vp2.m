% Calculate matrix of centrifugal and coriolis load on the joints for
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
% Cq [4x4]
%   matrix of coriolis and centrifugal joint torques

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:52
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function Cq = S4RRPP2_coriolismatJ_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP2_coriolismatJ_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From coriolismat_joint_fixb_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:52:25
% EndTime: 2018-11-14 13:52:25
% DurationCPUTime: 0.12s
% Computational Cost: add. (126->34), mult. (272->41), div. (0->0), fcn. (82->2), ass. (0->25)
t24 = mrSges(4,3) + mrSges(5,2);
t9 = (m(4) + m(5)) * qJ(3) + t24;
t37 = t9 * qJD(2);
t36 = t9 * qJD(3);
t35 = m(4) / 0.4e1 + m(5) / 0.4e1;
t34 = mrSges(3,1) + mrSges(4,1) + mrSges(5,1);
t33 = m(4) / 0.2e1;
t31 = m(5) / 0.2e1;
t16 = sin(qJ(2));
t28 = t16 * pkin(1);
t17 = cos(qJ(2));
t27 = t17 * pkin(1);
t12 = qJ(3) + t28;
t26 = t12 * t17;
t25 = t17 * mrSges(3,2);
t23 = t24 * t27;
t20 = -pkin(2) - t27;
t1 = -t23 + t34 * t28 + (t25 - m(4) * (t20 * t16 + t26) - m(5) * (t26 + (-pkin(3) + t20) * t16)) * pkin(1);
t22 = t1 * qJD(1);
t6 = 0.4e1 * t35 * t12 + t24;
t21 = t6 * qJD(1);
t18 = -qJD(1) * t9 - t37;
t13 = qJ(3) * t27;
t3 = (t33 + t31) * t28 + t24 + t35 * (0.4e1 * qJ(3) + 0.2e1 * t28);
t2 = [-t1 * qJD(2) + t6 * qJD(3), t3 * qJD(3) - t22 + (t23 + (-t34 * t16 - t25) * pkin(1) + 0.2e1 * (-pkin(2) * t28 + t13) * t33 + 0.2e1 * (t13 + (-pkin(2) - pkin(3)) * t28) * t31) * qJD(2), t3 * qJD(2) + t21, 0; t22 + t36, t36, -t18, 0; -t21 - t37, t18, 0, 0; 0, 0, 0, 0;];
Cq  = t2;
