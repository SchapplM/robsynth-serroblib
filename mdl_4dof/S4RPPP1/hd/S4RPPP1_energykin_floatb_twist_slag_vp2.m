% Calculate kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox (ehem. IRT-Maple-Toolbox)
% Datum: 2018-11-14 13:46
% Revision: ea61b7cc8771fdd0208f11149c97a676b461e858
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert( isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2018-11-14 13:45:22
% EndTime: 2018-11-14 13:45:23
% DurationCPUTime: 0.44s
% Computational Cost: add. (801->108), mult. (1251->140), div. (0->0), fcn. (918->6), ass. (0->37)
t46 = pkin(2) + qJ(4);
t30 = V_base(5) * pkin(5) + V_base(1);
t31 = -V_base(4) * pkin(5) + V_base(2);
t36 = sin(qJ(1));
t37 = cos(qJ(1));
t22 = t37 * t30 + t36 * t31;
t25 = t36 * V_base(5) + t37 * V_base(4);
t45 = qJ(2) * t25;
t44 = cos(pkin(6));
t24 = -t36 * V_base(4) + t37 * V_base(5);
t32 = V_base(6) + qJD(1);
t34 = sin(pkin(4));
t35 = cos(pkin(4));
t41 = t24 * t35 + t32 * t34;
t13 = t41 * qJ(2) + t22;
t21 = -t30 * t36 + t37 * t31;
t16 = pkin(1) * t32 - t35 * t45 + t21;
t19 = -pkin(1) * t24 - t34 * t45 + V_base(3);
t33 = sin(pkin(6));
t8 = t44 * t13 + (t16 * t35 + t19 * t34) * t33;
t43 = t34 * t44;
t42 = t35 * t44;
t9 = -t16 * t34 + t35 * t19 + qJD(2);
t20 = -t24 * t34 + t32 * t35;
t6 = -qJ(3) * t20 - t8;
t15 = t44 * t25 + t41 * t33;
t40 = -qJ(3) * t15 + t9;
t7 = -t33 * t13 + t16 * t42 + t19 * t43;
t39 = qJD(3) - t7;
t38 = V_base(3) ^ 2;
t14 = -t24 * t42 + t25 * t33 - t32 * t43;
t5 = -t20 * pkin(2) + t39;
t4 = pkin(2) * t14 + t40;
t3 = -pkin(3) * t14 + qJD(4) - t6;
t2 = t46 * t14 + t40;
t1 = t15 * pkin(3) - t46 * t20 + t39;
t10 = m(2) * (t21 ^ 2 + t22 ^ 2 + t38) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t38) / 0.2e1 + m(4) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(3) * (t7 ^ 2 + t8 ^ 2 + t9 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t21 * mrSges(2,1) - t22 * mrSges(2,2) + Ifges(2,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t21 * mrSges(2,3) + Ifges(2,5) * t32 + Ifges(2,1) * t25 / 0.2e1) * t25 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t22 * mrSges(2,3) + Ifges(2,4) * t25 + Ifges(2,6) * t32 + Ifges(2,2) * t24 / 0.2e1) * t24 + (t7 * mrSges(3,1) - t8 * mrSges(3,2) + t5 * mrSges(4,2) + t3 * mrSges(5,2) - t6 * mrSges(4,3) - t1 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * t20) * t20 + (t5 * mrSges(4,1) + t1 * mrSges(5,1) + t9 * mrSges(3,2) - t2 * mrSges(5,2) - t7 * mrSges(3,3) - t4 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(4,2) / 0.2e1) * t15 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t20) * t15 + (t9 * mrSges(3,1) + t6 * mrSges(4,1) - t3 * mrSges(5,1) - t4 * mrSges(4,2) - t8 * mrSges(3,3) + t2 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t14 + (Ifges(5,4) + Ifges(4,5) - Ifges(3,6)) * t20 + (-Ifges(3,4) - Ifges(4,6) + Ifges(5,6)) * t15) * t14;
T  = t10;
