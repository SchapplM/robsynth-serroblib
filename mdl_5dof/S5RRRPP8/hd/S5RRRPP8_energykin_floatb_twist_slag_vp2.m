% Calculate kinetic energy for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
% m_mdh [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:26
% EndTime: 2019-12-31 21:07:26
% DurationCPUTime: 0.59s
% Computational Cost: add. (1079->125), mult. (1364->161), div. (0->0), fcn. (964->6), ass. (0->40)
t49 = cos(qJ(3));
t48 = pkin(3) + qJ(5);
t42 = sin(qJ(1));
t44 = cos(qJ(1));
t32 = -t42 * V_base(4) + t44 * V_base(5);
t33 = t42 * V_base(5) + t44 * V_base(4);
t20 = -pkin(1) * t32 - pkin(6) * t33 + V_base(3);
t37 = V_base(5) * pkin(5) + V_base(1);
t38 = -V_base(4) * pkin(5) + V_base(2);
t29 = t44 * t37 + t42 * t38;
t39 = V_base(6) + qJD(1);
t24 = pkin(6) * t39 + t29;
t41 = sin(qJ(2));
t43 = cos(qJ(2));
t16 = t41 * t20 + t43 * t24;
t31 = qJD(2) - t32;
t13 = pkin(7) * t31 + t16;
t28 = -t42 * t37 + t38 * t44;
t23 = -pkin(1) * t39 - t28;
t26 = -t33 * t41 + t39 * t43;
t27 = t33 * t43 + t39 * t41;
t14 = -pkin(2) * t26 - pkin(7) * t27 + t23;
t40 = sin(qJ(3));
t8 = t49 * t13 + t40 * t14;
t15 = t20 * t43 - t41 * t24;
t25 = qJD(3) - t26;
t5 = -qJ(4) * t25 - t8;
t7 = -t40 * t13 + t49 * t14;
t12 = -pkin(2) * t31 - t15;
t47 = qJD(4) - t7;
t18 = t49 * t27 + t40 * t31;
t46 = -qJ(4) * t18 + t12;
t45 = V_base(3) ^ 2;
t17 = t27 * t40 - t49 * t31;
t6 = pkin(3) * t17 + t46;
t4 = -t25 * pkin(3) + t47;
t3 = t48 * t17 + t46;
t2 = -pkin(4) * t17 + qJD(5) - t5;
t1 = t18 * pkin(4) - t48 * t25 + t47;
t9 = m(2) * (t28 ^ 2 + t29 ^ 2 + t45) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t45) / 0.2e1 + m(3) * (t15 ^ 2 + t16 ^ 2 + t23 ^ 2) / 0.2e1 + m(4) * (t12 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(5) * (t4 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t28 * mrSges(2,1) - t29 * mrSges(2,2) + Ifges(2,3) * t39 / 0.2e1) * t39 + (t15 * mrSges(3,1) - t16 * mrSges(3,2) + Ifges(3,3) * t31 / 0.2e1) * t31 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t28 * mrSges(2,3) + Ifges(2,5) * t39 + Ifges(2,1) * t33 / 0.2e1) * t33 + (t23 * mrSges(3,2) - t15 * mrSges(3,3) + Ifges(3,5) * t31 + Ifges(3,1) * t27 / 0.2e1) * t27 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t29 * mrSges(2,3) + Ifges(2,4) * t33 + Ifges(2,6) * t39 + Ifges(2,2) * t32 / 0.2e1) * t32 + (-t23 * mrSges(3,1) + t16 * mrSges(3,3) + Ifges(3,4) * t27 + Ifges(3,6) * t31 + Ifges(3,2) * t26 / 0.2e1) * t26 + (t7 * mrSges(4,1) - t8 * mrSges(4,2) + t4 * mrSges(5,2) + t2 * mrSges(6,2) - t5 * mrSges(5,3) - t1 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t25) * t25 + (t4 * mrSges(5,1) + t1 * mrSges(6,1) + t12 * mrSges(4,2) - t3 * mrSges(6,2) - t7 * mrSges(4,3) - t6 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t18 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t25) * t18 + (t12 * mrSges(4,1) + t5 * mrSges(5,1) - t2 * mrSges(6,1) - t6 * mrSges(5,2) - t8 * mrSges(4,3) + t3 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t17 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t25 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t18) * t17;
T = t9;
