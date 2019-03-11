% Calculate kinetic energy for
% S6RPRRPP7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4]';
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
% Datum: 2019-03-09 04:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP7_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:50:11
% EndTime: 2019-03-09 04:50:12
% DurationCPUTime: 0.82s
% Computational Cost: add. (1327->146), mult. (1645->179), div. (0->0), fcn. (1120->6), ass. (0->46)
t61 = pkin(1) + pkin(7);
t60 = -pkin(4) - pkin(5);
t59 = cos(qJ(1));
t58 = cos(qJ(3));
t57 = cos(qJ(4));
t50 = sin(qJ(1));
t39 = t50 * V_base(5) + t59 * V_base(4);
t47 = V_base(6) + qJD(1);
t43 = V_base(5) * pkin(6) + V_base(1);
t44 = -V_base(4) * pkin(6) + V_base(2);
t33 = -t50 * t43 + t44 * t59;
t53 = qJD(2) - t33;
t19 = t39 * pkin(2) - t47 * t61 + t53;
t38 = t50 * V_base(4) - t59 * V_base(5);
t55 = -qJ(2) * t39 + V_base(3);
t24 = t38 * t61 + t55;
t49 = sin(qJ(3));
t14 = t49 * t19 + t58 * t24;
t37 = qJD(3) + t39;
t12 = pkin(8) * t37 + t14;
t34 = t59 * t43 + t50 * t44;
t29 = -t47 * qJ(2) - t34;
t25 = -pkin(2) * t38 - t29;
t31 = t58 * t38 - t47 * t49;
t32 = t49 * t38 + t47 * t58;
t16 = -pkin(3) * t31 - pkin(8) * t32 + t25;
t48 = sin(qJ(4));
t7 = t57 * t12 + t48 * t16;
t13 = t58 * t19 - t49 * t24;
t30 = qJD(4) - t31;
t5 = t30 * qJ(5) + t7;
t56 = pkin(3) * t37 + t13;
t6 = -t48 * t12 + t16 * t57;
t54 = qJD(5) - t6;
t23 = t32 * t57 + t48 * t37;
t52 = qJ(5) * t23 + t56;
t51 = V_base(3) ^ 2;
t27 = -t47 * pkin(1) + t53;
t26 = pkin(1) * t38 + t55;
t22 = t32 * t48 - t37 * t57;
t8 = pkin(4) * t22 - t52;
t4 = -t30 * pkin(4) + t54;
t3 = t22 * t60 + qJD(6) + t52;
t2 = qJ(6) * t22 + t5;
t1 = -t23 * qJ(6) + t60 * t30 + t54;
t9 = m(2) * (t33 ^ 2 + t34 ^ 2 + t51) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t51) / 0.2e1 + m(4) * (t13 ^ 2 + t14 ^ 2 + t25 ^ 2) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t29 ^ 2) / 0.2e1 + m(5) * (t56 ^ 2 + t6 ^ 2 + t7 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + m(6) * (t4 ^ 2 + t5 ^ 2 + t8 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t13 * mrSges(4,1) - t14 * mrSges(4,2) + Ifges(4,3) * t37 / 0.2e1) * t37 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t25 * mrSges(4,2) - t13 * mrSges(4,3) + Ifges(4,5) * t37 + Ifges(4,1) * t32 / 0.2e1) * t32 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t25 * mrSges(4,1) + t14 * mrSges(4,3) + Ifges(4,4) * t32 + Ifges(4,6) * t37 + Ifges(4,2) * t31 / 0.2e1) * t31 + (t33 * mrSges(2,1) - t34 * mrSges(2,2) + t27 * mrSges(3,2) - t29 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t47) * t47 + (t27 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t33 * mrSges(2,3) - t26 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t39 + (-Ifges(3,4) + Ifges(2,5)) * t47) * t39 + (V_base(3) * mrSges(2,1) + t29 * mrSges(3,1) - t26 * mrSges(3,2) - t34 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t38 + (Ifges(3,5) - Ifges(2,6)) * t47 + (-Ifges(2,4) - Ifges(3,6)) * t39) * t38 + (t6 * mrSges(5,1) - t4 * mrSges(6,1) - t1 * mrSges(7,1) - t7 * mrSges(5,2) + t2 * mrSges(7,2) + t5 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t30) * t30 + (-t56 * mrSges(5,2) + t4 * mrSges(6,2) + t3 * mrSges(7,2) - t6 * mrSges(5,3) - t8 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t23 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t30) * t23 + (-t56 * mrSges(5,1) + t8 * mrSges(6,1) - t3 * mrSges(7,1) - t5 * mrSges(6,2) - t7 * mrSges(5,3) + t2 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t22 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t30 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t23) * t22;
T  = t9;
