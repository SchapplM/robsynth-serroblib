% Calculate kinetic energy for
% S6RRPRRR14V3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
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
% Datum: 2019-04-12 15:12
% Revision: b693519ea345eb34ae9622239e7f1167217e9d53 (2019-04-09)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRR14V3_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(6,1),zeros(1,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRR14V3_energykin_floatb_twist_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-04-12 15:03:17
% EndTime: 2019-04-12 15:03:18
% DurationCPUTime: 1.02s
% Computational Cost: add. (1119->132), mult. (1515->189), div. (0->0), fcn. (1268->10), ass. (0->44)
t43 = sin(qJ(1));
t47 = cos(qJ(1));
t31 = t43 * V_base(5) + t47 * V_base(4);
t38 = V_base(6) + qJD(1);
t42 = sin(qJ(2));
t49 = cos(qJ(2));
t23 = t49 * t31 + t42 * t38;
t32 = t43 * V_base(1) - t47 * V_base(2);
t18 = -qJ(3) * t23 + t32;
t34 = t43 * V_base(2) + t47 * V_base(1);
t27 = t49 * t34 + t42 * V_base(3);
t30 = -t43 * V_base(4) + t47 * V_base(5);
t29 = qJD(2) - t30;
t19 = qJ(3) * t29 + t27;
t41 = sin(qJ(4));
t46 = cos(qJ(4));
t10 = t18 * t41 + t19 * t46;
t26 = -t42 * t34 + t49 * V_base(3);
t25 = qJD(3) - t26;
t40 = sin(qJ(5));
t45 = cos(qJ(5));
t5 = t10 * t40 - t45 * t25;
t52 = t5 ^ 2;
t8 = -t46 * t18 + t19 * t41;
t51 = t8 ^ 2;
t50 = t32 ^ 2;
t22 = t31 * t42 - t49 * t38;
t16 = t23 * t46 + t29 * t41;
t20 = qJD(4) + t22;
t12 = -t16 * t40 + t20 * t45;
t15 = -t23 * t41 + t29 * t46;
t48 = V_base(3) ^ 2;
t44 = cos(qJ(6));
t39 = sin(qJ(6));
t24 = t25 ^ 2;
t14 = qJD(5) - t15;
t13 = t16 * t45 + t20 * t40;
t11 = qJD(6) - t12;
t7 = t10 * t45 + t25 * t40;
t4 = t13 * t44 + t14 * t39;
t3 = -t13 * t39 + t14 * t44;
t2 = t39 * t8 + t44 * t7;
t1 = -t39 * t7 + t44 * t8;
t6 = (-t32 * mrSges(2,1) - t34 * mrSges(2,2) + Ifges(2,3) * t38 / 0.2e1) * t38 + (t25 * mrSges(5,2) + t8 * mrSges(5,3) + Ifges(5,5) * t20 + Ifges(5,1) * t16 / 0.2e1) * t16 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t26 * mrSges(3,1) - t25 * mrSges(4,1) - t27 * mrSges(3,2) + t19 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t29) * t29 + (V_base(3) * mrSges(2,2) + t32 * mrSges(2,3) + Ifges(2,5) * t38 + Ifges(2,1) * t31 / 0.2e1) * t31 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t8 * mrSges(6,2) + t5 * mrSges(6,3) + Ifges(6,5) * t14 + Ifges(6,1) * t13 / 0.2e1) * t13 + (-V_base(3) * mrSges(2,1) + t34 * mrSges(2,3) + Ifges(2,4) * t31 + Ifges(2,6) * t38 + Ifges(2,2) * t30 / 0.2e1) * t30 + (t5 * mrSges(7,2) - t1 * mrSges(7,3) + Ifges(7,1) * t4 / 0.2e1) * t4 + (-t25 * mrSges(5,1) + t10 * mrSges(5,3) + Ifges(5,4) * t16 + Ifges(5,6) * t20 + Ifges(5,2) * t15 / 0.2e1) * t15 + (-t8 * mrSges(6,1) + t7 * mrSges(6,3) + Ifges(6,4) * t13 + Ifges(6,6) * t14 + Ifges(6,2) * t12 / 0.2e1) * t12 + (t1 * mrSges(7,1) - t2 * mrSges(7,2) + Ifges(7,5) * t4 + Ifges(7,6) * t3 + Ifges(7,3) * t11 / 0.2e1) * t11 + (-t8 * mrSges(5,1) - t10 * mrSges(5,2) + Ifges(5,3) * t20 / 0.2e1) * t20 + m(2) * (t34 ^ 2 + t48 + t50) / 0.2e1 + m(3) * (t26 ^ 2 + t27 ^ 2 + t50) / 0.2e1 + m(5) * (t10 ^ 2 + t24 + t51) / 0.2e1 + m(6) * (t7 ^ 2 + t51 + t52) / 0.2e1 + m(7) * (t1 ^ 2 + t2 ^ 2 + t52) / 0.2e1 + (-t5 * mrSges(7,1) + t2 * mrSges(7,3) + Ifges(7,4) * t4 + Ifges(7,2) * t3 / 0.2e1) * t3 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t48) / 0.2e1 + m(4) * (t18 ^ 2 + t19 ^ 2 + t24) / 0.2e1 + (-t5 * mrSges(6,1) - t7 * mrSges(6,2) + Ifges(6,3) * t14 / 0.2e1) * t14 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t32 * mrSges(3,2) + t25 * mrSges(4,2) - t26 * mrSges(3,3) - t18 * mrSges(4,3) + (Ifges(4,1) / 0.2e1 + Ifges(3,1) / 0.2e1) * t23 + (Ifges(4,4) + Ifges(3,5)) * t29) * t23 + (t32 * mrSges(3,1) + t18 * mrSges(4,1) - t19 * mrSges(4,2) - t27 * mrSges(3,3) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t22 + (-Ifges(3,6) + Ifges(4,6)) * t29 + (-Ifges(3,4) + Ifges(4,5)) * t23) * t22;
T  = t6;
