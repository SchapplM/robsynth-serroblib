% Calculate kinetic energy for
% S6RPRPRR12
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
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6]';
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
% Datum: 2019-03-09 04:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR12_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR12_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:18:05
% EndTime: 2019-03-09 04:18:06
% DurationCPUTime: 0.92s
% Computational Cost: add. (1541->150), mult. (1905->196), div. (0->0), fcn. (1320->8), ass. (0->53)
t65 = pkin(1) + pkin(7);
t64 = pkin(3) + pkin(8);
t54 = sin(qJ(1));
t63 = cos(qJ(1));
t41 = t54 * V_base(4) - t63 * V_base(5);
t50 = V_base(6) + qJD(1);
t53 = sin(qJ(3));
t57 = cos(qJ(3));
t35 = t41 * t53 + t50 * t57;
t42 = t54 * V_base(5) + t63 * V_base(4);
t40 = qJD(3) + t42;
t46 = V_base(5) * pkin(6) + V_base(1);
t47 = -V_base(4) * pkin(6) + V_base(2);
t36 = -t54 * t46 + t47 * t63;
t60 = qJD(2) - t36;
t22 = t42 * pkin(2) - t50 * t65 + t60;
t62 = -qJ(2) * t42 + V_base(3);
t27 = t41 * t65 + t62;
t16 = t22 * t57 - t53 * t27;
t61 = qJD(4) - t16;
t10 = pkin(4) * t35 - t40 * t64 + t61;
t34 = -t57 * t41 + t50 * t53;
t37 = t63 * t46 + t54 * t47;
t31 = -t50 * qJ(2) - t37;
t28 = -pkin(2) * t41 - t31;
t59 = -qJ(4) * t35 + t28;
t13 = t34 * t64 + t59;
t52 = sin(qJ(5));
t56 = cos(qJ(5));
t6 = t52 * t10 + t56 * t13;
t17 = t53 * t22 + t57 * t27;
t15 = -t40 * qJ(4) - t17;
t5 = t56 * t10 - t13 * t52;
t11 = -pkin(4) * t34 - t15;
t33 = qJD(5) + t35;
t58 = V_base(3) ^ 2;
t55 = cos(qJ(6));
t51 = sin(qJ(6));
t32 = qJD(6) + t33;
t30 = -t50 * pkin(1) + t60;
t29 = pkin(1) * t41 + t62;
t26 = t34 * t52 + t40 * t56;
t25 = t34 * t56 - t40 * t52;
t20 = pkin(3) * t34 + t59;
t19 = t25 * t51 + t26 * t55;
t18 = t25 * t55 - t26 * t51;
t14 = -pkin(3) * t40 + t61;
t7 = -pkin(5) * t25 + t11;
t4 = pkin(9) * t25 + t6;
t3 = pkin(5) * t33 - pkin(9) * t26 + t5;
t2 = t3 * t51 + t4 * t55;
t1 = t3 * t55 - t4 * t51;
t8 = m(2) * (t36 ^ 2 + t37 ^ 2 + t58) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t58) / 0.2e1 + m(4) * (t16 ^ 2 + t17 ^ 2 + t28 ^ 2) / 0.2e1 + m(3) * (t29 ^ 2 + t30 ^ 2 + t31 ^ 2) / 0.2e1 + m(5) * (t14 ^ 2 + t15 ^ 2 + t20 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t5 * mrSges(6,1) - t6 * mrSges(6,2) + Ifges(6,3) * t33 / 0.2e1) * t33 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,3) * t32 / 0.2e1) * t32 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t11 * mrSges(6,2) - t5 * mrSges(6,3) + Ifges(6,5) * t33 + Ifges(6,1) * t26 / 0.2e1) * t26 + (t7 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,5) * t32 + Ifges(7,1) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-t11 * mrSges(6,1) + t6 * mrSges(6,3) + Ifges(6,4) * t26 + Ifges(6,6) * t33 + Ifges(6,2) * t25 / 0.2e1) * t25 + (-t7 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t19 + Ifges(7,6) * t32 + Ifges(7,2) * t18 / 0.2e1) * t18 + (t36 * mrSges(2,1) - t37 * mrSges(2,2) + t30 * mrSges(3,2) - t31 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t50) * t50 + (t16 * mrSges(4,1) - t17 * mrSges(4,2) + t14 * mrSges(5,2) - t15 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t40) * t40 + (t30 * mrSges(3,1) + V_base(3) * mrSges(2,2) - t36 * mrSges(2,3) - t29 * mrSges(3,3) + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1) * t42 + (-Ifges(3,4) + Ifges(2,5)) * t50) * t42 + (t14 * mrSges(5,1) + t28 * mrSges(4,2) - t16 * mrSges(4,3) - t20 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t35 + (-Ifges(5,4) + Ifges(4,5)) * t40) * t35 + (V_base(3) * mrSges(2,1) + t31 * mrSges(3,1) - t29 * mrSges(3,2) - t37 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t41 + (Ifges(3,5) - Ifges(2,6)) * t50 + (-Ifges(2,4) - Ifges(3,6)) * t42) * t41 + (t28 * mrSges(4,1) + t15 * mrSges(5,1) - t20 * mrSges(5,2) - t17 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t34 + (Ifges(5,5) - Ifges(4,6)) * t40 + (-Ifges(4,4) - Ifges(5,6)) * t35) * t34;
T  = t8;
