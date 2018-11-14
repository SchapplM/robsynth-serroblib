% Calculate kinetic energy for
% S3PRP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [3x1]
%   Generalized joint coordinates (joint angles)
% qJD [3x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [3x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,d2]';
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
% Datum: 2018-11-14 10:07
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S3PRP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(3,1),zeros(3,1),zeros(6,1),zeros(3,1),zeros(4,1),zeros(4,3),zeros(4,6)}
assert(isreal(qJ) && all(size(qJ) == [3 1]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: qJ has to be [3x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [3 1]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: qJD has to be [3x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [3 1]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: pkin has to be [3x1] (double)');
assert( isreal(m) && all(size(m) == [4 1]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: m has to be [4x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [4,3]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: mrSges has to be [4x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [4 6]), ...
  'S3PRP2_energykin_floatb_twist_slag_vp2: Ifges has to be [4x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 10:06:56
% EndTime: 2018-11-14 10:06:56
% DurationCPUTime: 0.24s
% Computational Cost: add. (187->76), mult. (240->92), div. (0->0), fcn. (72->2), ass. (0->17)
t19 = V_base(1) + qJD(1);
t10 = -V_base(6) * pkin(1) - V_base(5) * pkin(3) + t19;
t13 = V_base(6) * qJ(1) + V_base(2);
t11 = V_base(4) * pkin(3) + t13;
t20 = sin(qJ(2));
t21 = cos(qJ(2));
t5 = t20 * t10 + t21 * t11;
t14 = -V_base(5) * qJ(1) + V_base(3);
t4 = t21 * t10 - t20 * t11;
t12 = -V_base(4) * pkin(1) - t14;
t16 = -V_base(6) + qJD(2);
t9 = t20 * V_base(4) + t21 * V_base(5);
t8 = t20 * V_base(5) - t21 * V_base(4);
t3 = t16 * qJ(3) + t5;
t2 = -t16 * pkin(2) + qJD(3) - t4;
t1 = t8 * pkin(2) - t9 * qJ(3) + t12;
t6 = m(3) * (t12 ^ 2 + t4 ^ 2 + t5 ^ 2) / 0.2e1 + m(2) * (t13 ^ 2 + t14 ^ 2 + t19 ^ 2) / 0.2e1 + m(4) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - t19 * mrSges(2,1) - V_base(1) * mrSges(1,2) + t13 * mrSges(2,3) + (Ifges(1,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(6)) * V_base(6) + (t12 * mrSges(3,2) + t2 * mrSges(4,2) - t4 * mrSges(3,3) - t1 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t9) * t9 + (-V_base(3) * mrSges(1,1) + t19 * mrSges(2,2) + V_base(1) * mrSges(1,3) - t14 * mrSges(2,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(5) + (Ifges(2,4) + Ifges(1,6)) * V_base(6)) * V_base(5) + (t12 * mrSges(3,1) + t1 * mrSges(4,1) - t3 * mrSges(4,2) - t5 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t8 + (-Ifges(3,4) + Ifges(4,5)) * t9) * t8 + (t14 * mrSges(2,1) + V_base(3) * mrSges(1,2) - t13 * mrSges(2,2) - V_base(2) * mrSges(1,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(4) + (Ifges(1,5) + Ifges(2,6)) * V_base(6) + (Ifges(1,4) + Ifges(2,5)) * V_base(5)) * V_base(4) + (t4 * mrSges(3,1) - t2 * mrSges(4,1) - t5 * mrSges(3,2) + t3 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t16 + (Ifges(4,4) + Ifges(3,5)) * t9 + (-Ifges(3,6) + Ifges(4,6)) * t8) * t16;
T  = t6;
