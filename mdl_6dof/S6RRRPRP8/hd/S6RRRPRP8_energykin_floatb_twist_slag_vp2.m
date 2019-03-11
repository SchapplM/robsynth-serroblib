% Calculate kinetic energy for
% S6RRRPRP8
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
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP8_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP8_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:14:09
% EndTime: 2019-03-09 17:14:10
% DurationCPUTime: 0.92s
% Computational Cost: add. (2049->149), mult. (2573->196), div. (0->0), fcn. (1956->8), ass. (0->48)
t60 = sin(qJ(1));
t63 = cos(qJ(1));
t48 = -t60 * V_base(4) + t63 * V_base(5);
t49 = t60 * V_base(5) + t63 * V_base(4);
t33 = -pkin(1) * t48 - pkin(7) * t49 + V_base(3);
t54 = V_base(5) * pkin(6) + V_base(1);
t55 = -V_base(4) * pkin(6) + V_base(2);
t44 = t63 * t54 + t60 * t55;
t56 = V_base(6) + qJD(1);
t39 = pkin(7) * t56 + t44;
t59 = sin(qJ(2));
t62 = cos(qJ(2));
t26 = t59 * t33 + t62 * t39;
t47 = qJD(2) - t48;
t23 = pkin(8) * t47 + t26;
t43 = -t60 * t54 + t55 * t63;
t38 = -pkin(1) * t56 - t43;
t41 = -t49 * t59 + t62 * t56;
t42 = t49 * t62 + t56 * t59;
t24 = -pkin(2) * t41 - pkin(8) * t42 + t38;
t58 = sin(qJ(3));
t66 = cos(qJ(3));
t16 = t66 * t23 + t58 * t24;
t40 = qJD(3) - t41;
t13 = t40 * qJ(4) + t16;
t29 = t42 * t58 - t47 * t66;
t10 = pkin(9) * t29 + t13;
t57 = sin(qJ(5));
t61 = cos(qJ(5));
t30 = t42 * t66 + t58 * t47;
t15 = -t58 * t23 + t24 * t66;
t65 = qJD(4) - t15;
t8 = -t30 * pkin(9) + (-pkin(3) - pkin(4)) * t40 + t65;
t4 = t61 * t10 + t57 * t8;
t25 = t62 * t33 - t59 * t39;
t22 = -t47 * pkin(2) - t25;
t3 = -t10 * t57 + t61 * t8;
t14 = t29 * pkin(3) - t30 * qJ(4) + t22;
t11 = -pkin(4) * t29 - t14;
t64 = V_base(3) ^ 2;
t37 = qJD(5) - t40;
t18 = t29 * t57 + t30 * t61;
t17 = t29 * t61 - t30 * t57;
t12 = -t40 * pkin(3) + t65;
t5 = -pkin(5) * t17 + qJD(6) + t11;
t2 = qJ(6) * t17 + t4;
t1 = pkin(5) * t37 - qJ(6) * t18 + t3;
t6 = m(2) * (t43 ^ 2 + t44 ^ 2 + t64) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t64) / 0.2e1 + m(3) * (t25 ^ 2 + t26 ^ 2 + t38 ^ 2) / 0.2e1 + m(4) * (t15 ^ 2 + t16 ^ 2 + t22 ^ 2) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(6) * (t11 ^ 2 + t3 ^ 2 + t4 ^ 2) / 0.2e1 + m(5) * (t12 ^ 2 + t13 ^ 2 + t14 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t43 * mrSges(2,1) - t44 * mrSges(2,2) + Ifges(2,3) * t56 / 0.2e1) * t56 + (t25 * mrSges(3,1) - t26 * mrSges(3,2) + Ifges(3,3) * t47 / 0.2e1) * t47 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (V_base(3) * mrSges(2,2) - t43 * mrSges(2,3) + Ifges(2,5) * t56 + Ifges(2,1) * t49 / 0.2e1) * t49 + (t38 * mrSges(3,2) - t25 * mrSges(3,3) + Ifges(3,5) * t47 + Ifges(3,1) * t42 / 0.2e1) * t42 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (-V_base(3) * mrSges(2,1) + t44 * mrSges(2,3) + Ifges(2,4) * t49 + Ifges(2,6) * t56 + Ifges(2,2) * t48 / 0.2e1) * t48 + (-t38 * mrSges(3,1) + t26 * mrSges(3,3) + Ifges(3,4) * t42 + Ifges(3,6) * t47 + Ifges(3,2) * t41 / 0.2e1) * t41 + (t15 * mrSges(4,1) - t12 * mrSges(5,1) - t16 * mrSges(4,2) + t13 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t40) * t40 + (t3 * mrSges(6,1) + t1 * mrSges(7,1) - t4 * mrSges(6,2) - t2 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t37) * t37 + (t22 * mrSges(4,2) + t12 * mrSges(5,2) - t15 * mrSges(4,3) - t14 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t30 + (Ifges(5,4) + Ifges(4,5)) * t40) * t30 + (t11 * mrSges(6,2) + t5 * mrSges(7,2) - t3 * mrSges(6,3) - t1 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t18 + (Ifges(6,5) + Ifges(7,5)) * t37) * t18 + (t22 * mrSges(4,1) + t14 * mrSges(5,1) - t13 * mrSges(5,2) - t16 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t29 + (-Ifges(4,6) + Ifges(5,6)) * t40 + (-Ifges(4,4) + Ifges(5,5)) * t30) * t29 + (-t11 * mrSges(6,1) - t5 * mrSges(7,1) + t4 * mrSges(6,3) + t2 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t17 + (Ifges(6,6) + Ifges(7,6)) * t37 + (Ifges(6,4) + Ifges(7,4)) * t18) * t17;
T  = t6;
