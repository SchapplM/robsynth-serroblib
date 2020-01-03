% Calculate kinetic energy for
% S4RPPR7
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
%   pkin=[a2,a3,a4,d1,d4,theta3]';
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
% Datum: 2019-12-31 16:41
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPR7_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPR7_energykin_floatb_twist_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:41:33
% EndTime: 2019-12-31 16:41:34
% DurationCPUTime: 0.44s
% Computational Cost: add. (593->105), mult. (775->141), div. (0->0), fcn. (480->6), ass. (0->36)
t43 = cos(qJ(1));
t42 = pkin(1) + qJ(3);
t37 = sin(qJ(1));
t25 = t37 * V_base(5) + t43 * V_base(4);
t33 = V_base(6) + qJD(1);
t29 = V_base(5) * pkin(4) + V_base(1);
t30 = -V_base(4) * pkin(4) + V_base(2);
t21 = -t37 * t29 + t43 * t30;
t40 = qJD(2) - t21;
t12 = t25 * pkin(2) - t42 * t33 + t40;
t24 = t37 * V_base(4) - t43 * V_base(5);
t41 = -qJ(2) * t25 + V_base(3);
t14 = t42 * t24 + t41;
t34 = sin(pkin(6));
t35 = cos(pkin(6));
t6 = t34 * t12 + t35 * t14;
t22 = t43 * t29 + t37 * t30;
t18 = -t33 * qJ(2) - t22;
t5 = t35 * t12 - t14 * t34;
t15 = -pkin(2) * t24 + qJD(3) - t18;
t39 = V_base(3) ^ 2;
t38 = cos(qJ(4));
t36 = sin(qJ(4));
t23 = qJD(4) + t25;
t20 = t24 * t34 + t33 * t35;
t19 = t24 * t35 - t33 * t34;
t17 = -t33 * pkin(1) + t40;
t16 = pkin(1) * t24 + t41;
t9 = t19 * t36 + t20 * t38;
t8 = t19 * t38 - t20 * t36;
t7 = -pkin(3) * t19 + t15;
t4 = pkin(5) * t19 + t6;
t3 = pkin(3) * t25 - pkin(5) * t20 + t5;
t2 = t3 * t36 + t38 * t4;
t1 = t3 * t38 - t36 * t4;
t10 = m(2) * (t21 ^ 2 + t22 ^ 2 + t39) / 0.2e1 + m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + t39) / 0.2e1 + m(4) * (t15 ^ 2 + t5 ^ 2 + t6 ^ 2) / 0.2e1 + m(3) * (t16 ^ 2 + t17 ^ 2 + t18 ^ 2) / 0.2e1 + m(5) * (t1 ^ 2 + t2 ^ 2 + t7 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) - V_base(1) * mrSges(1,2) + Ifges(1,3) * V_base(6) / 0.2e1) * V_base(6) + (t7 * mrSges(5,2) - t1 * mrSges(5,3) + Ifges(5,1) * t9 / 0.2e1) * t9 + (t15 * mrSges(4,2) - t5 * mrSges(4,3) + Ifges(4,1) * t20 / 0.2e1) * t20 + (-V_base(3) * mrSges(1,1) + V_base(1) * mrSges(1,3) + Ifges(1,6) * V_base(6) + Ifges(1,2) * V_base(5) / 0.2e1) * V_base(5) + (-t7 * mrSges(5,1) + t2 * mrSges(5,3) + Ifges(5,4) * t9 + Ifges(5,2) * t8 / 0.2e1) * t8 + (-t15 * mrSges(4,1) + t6 * mrSges(4,3) + Ifges(4,4) * t20 + Ifges(4,2) * t19 / 0.2e1) * t19 + (V_base(3) * mrSges(1,2) - V_base(2) * mrSges(1,3) + Ifges(1,4) * V_base(5) + Ifges(1,5) * V_base(6) + Ifges(1,1) * V_base(4) / 0.2e1) * V_base(4) + (t1 * mrSges(5,1) - t2 * mrSges(5,2) + Ifges(5,5) * t9 + Ifges(5,6) * t8 + Ifges(5,3) * t23 / 0.2e1) * t23 + (t21 * mrSges(2,1) - t22 * mrSges(2,2) + t17 * mrSges(3,2) - t18 * mrSges(3,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t33) * t33 + (V_base(3) * mrSges(2,1) + t18 * mrSges(3,1) - t16 * mrSges(3,2) - t22 * mrSges(2,3) + (Ifges(2,2) / 0.2e1 + Ifges(3,3) / 0.2e1) * t24 + (Ifges(3,5) - Ifges(2,6)) * t33) * t24 + (t17 * mrSges(3,1) + t5 * mrSges(4,1) + V_base(3) * mrSges(2,2) - t6 * mrSges(4,2) - t21 * mrSges(2,3) - t16 * mrSges(3,3) + Ifges(4,5) * t20 + Ifges(4,6) * t19 + (Ifges(2,1) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t25 + (-Ifges(3,4) + Ifges(2,5)) * t33 + (-Ifges(2,4) - Ifges(3,6)) * t24) * t25;
T = t10;
