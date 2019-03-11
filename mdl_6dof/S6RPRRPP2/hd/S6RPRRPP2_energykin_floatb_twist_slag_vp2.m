% Calculate kinetic energy for
% S6RPRRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP2_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP2_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:31:28
% EndTime: 2019-03-09 04:31:29
% DurationCPUTime: 0.93s
% Computational Cost: add. (2087->149), mult. (3001->194), div. (0->0), fcn. (2340->8), ass. (0->49)
t65 = -pkin(4) - pkin(5);
t64 = cos(qJ(3));
t63 = cos(qJ(4));
t50 = V_base(5) * pkin(6) + V_base(1);
t51 = -V_base(4) * pkin(6) + V_base(2);
t57 = sin(qJ(1));
t58 = cos(qJ(1));
t41 = -t50 * t57 + t58 * t51;
t45 = t57 * V_base(5) + t58 * V_base(4);
t52 = V_base(6) + qJD(1);
t33 = pkin(1) * t52 - qJ(2) * t45 + t41;
t42 = t58 * t50 + t57 * t51;
t44 = -t57 * V_base(4) + t58 * V_base(5);
t37 = qJ(2) * t44 + t42;
t53 = sin(pkin(9));
t54 = cos(pkin(9));
t27 = t53 * t33 + t54 * t37;
t20 = pkin(7) * t52 + t27;
t39 = t44 * t54 - t45 * t53;
t40 = t44 * t53 + t45 * t54;
t43 = -pkin(1) * t44 + qJD(2) + V_base(3);
t23 = -pkin(2) * t39 - pkin(7) * t40 + t43;
t56 = sin(qJ(3));
t15 = t64 * t20 + t56 * t23;
t38 = qJD(3) - t39;
t12 = pkin(8) * t38 + t15;
t26 = t33 * t54 - t53 * t37;
t19 = -pkin(2) * t52 - t26;
t31 = -t40 * t56 + t64 * t52;
t32 = t64 * t40 + t56 * t52;
t16 = -pkin(3) * t31 - pkin(8) * t32 + t19;
t55 = sin(qJ(4));
t7 = t63 * t12 + t55 * t16;
t14 = -t56 * t20 + t64 * t23;
t29 = qJD(4) - t31;
t5 = t29 * qJ(5) + t7;
t62 = pkin(3) * t38 + t14;
t6 = -t55 * t12 + t63 * t16;
t61 = qJD(5) - t6;
t25 = t63 * t32 + t55 * t38;
t60 = qJ(5) * t25 + t62;
t59 = V_base(3) ^ 2;
t24 = t32 * t55 - t63 * t38;
t8 = pkin(4) * t24 - t60;
t4 = -t29 * pkin(4) + t61;
t3 = t65 * t24 + qJD(6) + t60;
t2 = qJ(6) * t24 + t5;
t1 = -t25 * qJ(6) + t65 * t29 + t61;
t9 = m(2) * (t41 ^ 2 + t42 ^ 2 + t59) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t59) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t43 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t62 ^ 2 + t7 ^ 2) / 0.2e1 + m(4) * (t14 ^ 2 + t15 ^ 2 + t19 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (V_base(3) * mrSges(2,2) - t41 * mrSges(2,3) + Ifges(2,1) * t45 / 0.2e1) * t45 + (t43 * mrSges(3,2) - t26 * mrSges(3,3) + Ifges(3,1) * t40 / 0.2e1) * t40 + (t14 * mrSges(4,1) - t15 * mrSges(4,2) + Ifges(4,3) * t38 / 0.2e1) * t38 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-V_base(3) * mrSges(2,1) + t42 * mrSges(2,3) + Ifges(2,4) * t45 + Ifges(2,2) * t44 / 0.2e1) * t44 + (-t43 * mrSges(3,1) + t27 * mrSges(3,3) + Ifges(3,4) * t40 + Ifges(3,2) * t39 / 0.2e1) * t39 + (t19 * mrSges(4,2) - t14 * mrSges(4,3) + Ifges(4,5) * t38 + Ifges(4,1) * t32 / 0.2e1) * t32 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t19 * mrSges(4,1) + t15 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t38 + Ifges(4,2) * t31 / 0.2e1) * t31 + (t41 * mrSges(2,1) + t26 * mrSges(3,1) - t42 * mrSges(2,2) - t27 * mrSges(3,2) + Ifges(2,5) * t45 + Ifges(3,5) * t40 + Ifges(2,6) * t44 + Ifges(3,6) * t39 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * t52) * t52 + (t6 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t5 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29) * t29 + (-t62 * mrSges(5,2) + t4 * mrSges(6,2) + t3 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t25 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t29) * t25 + (-t62 * mrSges(5,1) + t8 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) - t7 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t24 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t29 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t25) * t24;
T  = t9;
