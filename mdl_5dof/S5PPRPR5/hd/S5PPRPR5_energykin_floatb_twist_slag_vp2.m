% Calculate kinetic energy for
% S5PPRPR5
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
%   pkin=[a2,a3,a4,a5,d3,d5,theta1]';
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
% Datum: 2019-12-31 17:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRPR5_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRPR5_energykin_floatb_twist_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:33:22
% EndTime: 2019-12-31 17:33:23
% DurationCPUTime: 0.57s
% Computational Cost: add. (733->126), mult. (1060->160), div. (0->0), fcn. (676->6), ass. (0->40)
t50 = pkin(3) + pkin(6);
t49 = cos(qJ(3));
t41 = sin(pkin(7));
t48 = cos(pkin(7));
t31 = t41 * V_base(5) + t48 * V_base(4);
t35 = V_base(5) * qJ(1) + V_base(1);
t36 = -V_base(4) * qJ(1) + V_base(2);
t25 = -t41 * t35 + t36 * t48;
t46 = qJD(2) - t25;
t13 = -t31 * pkin(5) + (-pkin(1) - pkin(2)) * V_base(6) + t46;
t26 = t48 * t35 + t41 * t36;
t24 = V_base(6) * qJ(2) + t26;
t30 = t41 * V_base(4) - t48 * V_base(5);
t18 = pkin(5) * t30 + t24;
t43 = sin(qJ(3));
t10 = t43 * t13 + t49 * t18;
t40 = V_base(3) + qJD(1);
t38 = -V_base(6) + qJD(3);
t8 = -qJ(4) * t38 - t10;
t9 = t13 * t49 - t43 * t18;
t20 = t30 * pkin(1) - t31 * qJ(2) + t40;
t47 = qJD(4) - t9;
t22 = t43 * t30 + t31 * t49;
t12 = -pkin(2) * t30 - t20;
t45 = -qJ(4) * t22 + t12;
t44 = cos(qJ(5));
t42 = sin(qJ(5));
t23 = -V_base(6) * pkin(1) + t46;
t21 = -t30 * t49 + t31 * t43;
t19 = qJD(5) + t22;
t15 = t42 * t21 + t38 * t44;
t14 = t21 * t44 - t42 * t38;
t7 = -t38 * pkin(3) + t47;
t6 = pkin(3) * t21 + t45;
t5 = -pkin(4) * t21 - t8;
t4 = t22 * pkin(4) - t38 * t50 + t47;
t3 = t21 * t50 + t45;
t2 = t3 * t44 + t42 * t4;
t1 = -t42 * t3 + t4 * t44;
t11 = m(2) * (t25 ^ 2 + t26 ^ 2 + t40 ^ 2) / 0.2e1 + m(3) * (t20 ^ 2 + t23 ^ 2 + t24 ^ 2) / 0.2e1 + m(5) * (t6 ^ 2 + t7 ^ 2 + t8 ^ 2) / 0.2e1 + m(4) * (t10 ^ 2 + t12 ^ 2 + t9 ^ 2) / 0.2e1 + m(6) * (t1 ^ 2 + t2 ^ 2 + t5 ^ 2) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (t1 * mrSges(6,1) - t2 * mrSges(6,2) + Ifges(6,3) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t5 * mrSges(6,2) - t1 * mrSges(6,3) + Ifges(6,5) * t19 + Ifges(6,1) * t15 / 0.2e1) * t15 + (-t5 * mrSges(6,1) + t2 * mrSges(6,3) + Ifges(6,4) * t15 + Ifges(6,6) * t19 + Ifges(6,2) * t14 / 0.2e1) * t14 + (t9 * mrSges(4,1) - t10 * mrSges(4,2) + t7 * mrSges(5,2) - t8 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t38) * t38 + (t40 * mrSges(2,2) + t23 * mrSges(3,2) - t25 * mrSges(2,3) - t20 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * t31) * t31 + (t40 * mrSges(2,1) + t20 * mrSges(3,1) - t24 * mrSges(3,2) - t26 * mrSges(2,3) + (Ifges(3,3) / 0.2e1 + Ifges(2,2) / 0.2e1) * t30 + (-Ifges(2,4) + Ifges(3,5)) * t31) * t30 + (t7 * mrSges(5,1) + t12 * mrSges(4,2) - t9 * mrSges(4,3) - t6 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t22 + (-Ifges(5,4) + Ifges(4,5)) * t38) * t22 + (t12 * mrSges(4,1) + t8 * mrSges(5,1) - t6 * mrSges(5,2) - t10 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t21 + (Ifges(5,5) - Ifges(4,6)) * t38 + (-Ifges(4,4) - Ifges(5,6)) * t22) * t21 + (V_base(2) * mrSges(1,1) + t25 * mrSges(2,1) - t23 * mrSges(3,1) - V_base(1) * mrSges(1,2) - t26 * mrSges(2,2) + t24 * mrSges(3,3) + Ifges(1,5) * V_base(4) + Ifges(1,6) * V_base(5) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(1,3) / 0.2e1) * V_base(6) + (Ifges(3,4) + Ifges(2,5)) * t31 + (-Ifges(2,6) + Ifges(3,6)) * t30) * V_base(6);
T = t11;
