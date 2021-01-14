% Calculate kinetic energy for
% S1P1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [1x1]
%   Generalized joint coordinates (joint angles)
% qJD [1x1]
%   Generalized joint velocities
% V_base [6x1]
%   Base Velocity (twist: stacked translational and angular velocity) in base frame
% pkin [1x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[dummy]';
% m [2x1]
%   mass of all robot links (including the base)
% mrSges [2x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [2x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2021-01-14 12:22
% Revision: 96facaeb42edba38506bd76ea342a8981e82f256 (2020-11-21)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S1P1_energykin_floatb_twist_slag_vp2(qJ, qJD, V_base, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(1,1),zeros(1,1),zeros(6,1),zeros(1,1),zeros(2,1),zeros(2,3),zeros(2,6)}
assert(isreal(qJ) && all(size(qJ) == [1 1]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: qJ has to be [1x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [1 1]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: qJD has to be [1x1] (double)');
assert(isreal(V_base) && all(size(V_base) == [6 1]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: V_base has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [1 1]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: pkin has to be [1x1] (double)');
assert(isreal(m) && all(size(m) == [2 1]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: m has to be [2x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [2,3]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: mrSges has to be [2x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [2 6]), ...
  'S1P1_energykin_floatb_twist_slag_vp2: Ifges has to be [2x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_floatb_twist_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2021-01-14 12:21:50
% EndTime: 2021-01-14 12:21:50
% DurationCPUTime: 0.14s
% Computational Cost: add. (41->32), mult. (62->39), div. (0->0), fcn. (0->0), ass. (0->4)
t3 = V_base(3) + qJD(1);
t2 = -V_base(4) * qJ(1) + V_base(2);
t1 = V_base(5) * qJ(1) + V_base(1);
t4 = m(1) * (V_base(1) ^ 2 + V_base(2) ^ 2 + V_base(3) ^ 2) / 0.2e1 + m(2) * (t1 ^ 2 + t2 ^ 2 + t3 ^ 2) / 0.2e1 + (V_base(2) * mrSges(1,1) + t2 * mrSges(2,1) - V_base(1) * mrSges(1,2) - t1 * mrSges(2,2) + (Ifges(1,3) / 0.2e1 + Ifges(2,3) / 0.2e1) * V_base(6)) * V_base(6) + (-V_base(3) * mrSges(1,1) - t3 * mrSges(2,1) + V_base(1) * mrSges(1,3) + t1 * mrSges(2,3) + (Ifges(1,2) / 0.2e1 + Ifges(2,2) / 0.2e1) * V_base(5) + (Ifges(1,6) + Ifges(2,6)) * V_base(6)) * V_base(5) + (V_base(3) * mrSges(1,2) + t3 * mrSges(2,2) - V_base(2) * mrSges(1,3) - t2 * mrSges(2,3) + (Ifges(1,1) / 0.2e1 + Ifges(2,1) / 0.2e1) * V_base(4) + (Ifges(1,5) + Ifges(2,5)) * V_base(6) + (Ifges(1,4) + Ifges(2,4)) * V_base(5)) * V_base(4);
T = t4;
