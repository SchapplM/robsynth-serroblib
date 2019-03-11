% Calculate kinetic energy for
% S4PPRR2
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
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-08 18:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4PPRR2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4PPRR2_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:16:48
% EndTime: 2019-03-08 18:16:49
% DurationCPUTime: 0.40s
% Computational Cost: add. (563->104), mult. (795->142), div. (0->0), fcn. (468->6), ass. (0->33)
t33 = V_base(2) + qJD(1);
t25 = V_base(6) * pkin(1) - V_base(4) * qJ(2) + t33;
t28 = -V_base(6) * qJ(1) + V_base(1);
t27 = V_base(5) * qJ(2) + t28;
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t17 = t35 * t25 - t27 * t34;
t24 = t34 * V_base(5) + t35 * V_base(4);
t12 = V_base(6) * pkin(2) - pkin(4) * t24 + t17;
t18 = t34 * t25 + t35 * t27;
t23 = -t34 * V_base(4) + t35 * V_base(5);
t14 = pkin(4) * t23 + t18;
t37 = sin(qJ(3));
t39 = cos(qJ(3));
t6 = t37 * t12 + t39 * t14;
t29 = -V_base(4) * qJ(1) - V_base(3);
t5 = t39 * t12 - t14 * t37;
t31 = V_base(6) + qJD(3);
t26 = -V_base(5) * pkin(1) + qJD(2) - t29;
t19 = -pkin(2) * t23 + t26;
t38 = cos(qJ(4));
t36 = sin(qJ(4));
t30 = qJD(4) + t31;
t16 = t37 * t23 + t24 * t39;
t15 = t23 * t39 - t37 * t24;
t9 = -pkin(3) * t15 + t19;
t8 = t15 * t36 + t16 * t38;
t7 = t15 * t38 - t16 * t36;
t4 = pkin(5) * t15 + t6;
t3 = pkin(3) * t31 - pkin(5) * t16 + t5;
t2 = t3 * t36 + t38 * t4;
t1 = t3 * t38 - t36 * t4;
t10 = m(2) * (t28 ^ 2 + t29 ^ 2 + t33 ^ 2) / 0.2e1 + m(3) * (t17 ^ 2 + t18 ^ 2 + t26 ^ 2) / 0.2e1 + m(4) * (t19 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t9 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (t9 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t8 / 0.2e1) * t8 + (t5 * mrSges(4,1) - t6 * mrSges(4,2) + Ifges(4,3) * t31 / 0.2e1) * t31 + (t26 * mrSges(3,2) - t17 * mrSges(3,3) + Ifges(3,1) * t24 / 0.2e1) * t24 + (-t9 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t8 + Ifges(5,2) * t7 / 0.2e1) * t7 + (-t26 * mrSges(3,1) + t18 * mrSges(3,3) + Ifges(3,4) * t24 + Ifges(3,2) * t23 / 0.2e1) * t23 + (t19 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,5) * t31 + Ifges(4,1) * t16 / 0.2e1) * t16 + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t8 + Ifges(5,6) * t7 + Ifges(5,3) * t30 / 0.2e1) * t30 + (-t19 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t16 + Ifges(4,6) * t31 + Ifges(4,2) * t15 / 0.2e1) * t15 + (-V_base(3) * mrSges(1,1) + t29 * mrSges(2,1) - t28 * mrSges(2,2) + V_base(1) * mrSges(1,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(5)) * V_base(5) + (V_base(3) * mrSges(1,2) + t33 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t29 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(4) + (Ifges(1,4) + Ifges(2,5)) * V_base(5)) * V_base(4) + (V_base(2) * mrSges(1,1) + t33 * mrSges(2,1) + t17 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t18 * mrSges(3,2) - t28 * mrSges(2,3) + Ifges(3,5) * t24 + Ifges(3,6) * t23 + (Ifges(3,3) / 0.2e1 + Ifges(1,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(6) + (Ifges(1,6) - Ifges(2,6)) * V_base(5) + (-Ifges(2,4) + Ifges(1,5)) * V_base(4)) * V_base(6);
T  = t10;
