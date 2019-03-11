% Calculate kinetic energy for
% S4RPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,theta2]';
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
% Datum: 2019-03-08 18:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRP1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:29:26
% EndTime: 2019-03-08 18:29:26
% DurationCPUTime: 0.06s
% Computational Cost: add. (50->25), mult. (111->39), div. (0->0), fcn. (36->4), ass. (0->14)
t59 = m(3) / 0.2e1;
t52 = cos(pkin(6));
t47 = (pkin(1) * t52 + pkin(2)) * qJD(1);
t53 = sin(qJ(3));
t54 = cos(qJ(3));
t51 = sin(pkin(6));
t58 = pkin(1) * qJD(1) * t51;
t45 = t53 * t47 + t54 * t58;
t44 = t54 * t47 - t53 * t58;
t55 = qJD(2) ^ 2;
t50 = qJD(1) + qJD(3);
t43 = t50 * qJ(4) + t45;
t42 = -t50 * pkin(3) + qJD(4) - t44;
t1 = m(4) * (t44 ^ 2 + t45 ^ 2 + t55) / 0.2e1 + m(5) * (t42 ^ 2 + t43 ^ 2 + t55) / 0.2e1 + t55 * t59 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t52 * mrSges(3,1) - t51 * mrSges(3,2) + (t51 ^ 2 + t52 ^ 2) * t59 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (t44 * mrSges(4,1) - t42 * mrSges(5,1) - t45 * mrSges(4,2) + t43 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t50) * t50;
T  = t1;
