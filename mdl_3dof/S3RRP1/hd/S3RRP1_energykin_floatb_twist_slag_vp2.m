% Calculate kinetic energy for
% S3RRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [4x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d1,d2]';
% m_mdh [4x1]
%   mass of all robot links (including the base)
% mrSges [4x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [4x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 10:15
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3RRP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(4,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [4 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: pkin has to be [4x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3RRP1_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:15:06
% EndTime: 2018-11-14 10:15:06
% DurationCPUTime: 0.28s
% Computational Cost: add. (339->80), mult. (480->108), div. (0->0), fcn. (288->4), ass. (0->24)
t25 = sin(qJ(2));
t29 = cos(qJ(2));
t21 = V_base(5) * pkin(3) + V_base(1);
t22 = -V_base(4) * pkin(3) + V_base(2);
t26 = sin(qJ(1));
t27 = cos(qJ(1));
t12 = -t26 * t21 + t27 * t22;
t17 = t26 * V_base(5) + t27 * V_base(4);
t24 = V_base(6) + qJD(1);
t7 = t24 * pkin(1) - t17 * pkin(4) + t12;
t13 = t27 * t21 + t26 * t22;
t16 = -t26 * V_base(4) + t27 * V_base(5);
t9 = t16 * pkin(4) + t13;
t5 = t25 * t7 + t29 * t9;
t14 = -t16 * pkin(1) + V_base(3);
t4 = -t25 * t9 + t29 * t7;
t28 = V_base(3) ^ 2;
t23 = qJD(2) + t24;
t11 = t25 * t16 + t17 * t29;
t10 = -t16 * t29 + t25 * t17;
t3 = t10 * pkin(2) - t11 * qJ(3) + t14;
t2 = t23 * qJ(3) + t5;
t1 = -t23 * pkin(2) + qJD(3) - t4;
t6 = m(2) * (t12 ^ 2 + t13 ^ 2 + t28) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t28) / 0.2e1 + m(3) * (t14 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + m(4) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t12 * mrSges(2,1) - t13 * mrSges(2,2) + Ifges(2,3) * t24 / 0.2e1) * t24 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t12 * mrSges(2,3) + Ifges(2,5) * t24 + Ifges(2,1) * t17 / 0.2e1) * t17 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t13 * mrSges(2,3) + Ifges(2,4) * t17 + Ifges(2,6) * t24 + Ifges(2,2) * t16 / 0.2e1) * t16 + (t4 * mrSges(3,1) - t1 * mrSges(4,1) - t5 * mrSges(3,2) + t2 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t23) * t23 + (t14 * mrSges(3,2) + t1 * mrSges(4,2) - t4 * mrSges(3,3) - t3 * mrSges(4,3) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t11 + (Ifges(4,4) + Ifges(3,5)) * t23) * t11 + (t14 * mrSges(3,1) + t3 * mrSges(4,1) - t2 * mrSges(4,2) - t5 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t10 + (-Ifges(3,6) + Ifges(4,6)) * t23 + (-Ifges(3,4) + Ifges(4,5)) * t11) * t10;
T  = t6;
