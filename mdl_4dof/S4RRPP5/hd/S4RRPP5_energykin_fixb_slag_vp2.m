% Calculate kinetic energy for
% S4RRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d2]';
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
% Datum: 2019-12-31 17:00
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRPP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(5,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRPP5_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:00:07
% EndTime: 2019-12-31 17:00:08
% DurationCPUTime: 0.18s
% Computational Cost: add. (84->57), mult. (199->71), div. (0->0), fcn. (62->2), ass. (0->16)
t54 = pkin(5) * mrSges(3,3);
t53 = -pkin(2) - qJ(4);
t44 = sin(qJ(2));
t52 = qJD(1) * t44;
t45 = cos(qJ(2));
t51 = qJD(1) * t45;
t50 = pkin(5) * t52 + qJD(3);
t49 = qJD(2) * qJ(3);
t48 = -qJ(3) * t44 - pkin(1);
t42 = -pkin(5) * t51 - t49;
t41 = -qJD(2) * pkin(2) + t50;
t40 = (-pkin(2) * t45 + t48) * qJD(1);
t39 = t49 + qJD(4) + (pkin(3) + pkin(5)) * t51;
t38 = pkin(3) * t52 + t53 * qJD(2) + t50;
t37 = (t53 * t45 + t48) * qJD(1);
t1 = m(4) * (t40 ^ 2 + t41 ^ 2 + t42 ^ 2) / 0.2e1 + m(5) * (t37 ^ 2 + t38 ^ 2 + t39 ^ 2) / 0.2e1 + (t41 * mrSges(4,2) + t39 * mrSges(5,2) - t42 * mrSges(4,3) - t38 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + ((-t42 * mrSges(4,1) + t39 * mrSges(5,1) + t40 * mrSges(4,2) - t37 * mrSges(5,3) + (-mrSges(3,2) * pkin(5) - Ifges(5,4) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t45 + (t41 * mrSges(4,1) + t38 * mrSges(5,1) - t37 * mrSges(5,2) - t40 * mrSges(4,3) + (-mrSges(3,1) * pkin(5) - Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * qJD(2)) * t44 + (m(3) * (pkin(1) ^ 2 + (t44 ^ 2 + t45 ^ 2) * pkin(5) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t54 + Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t45) * t45 + (-pkin(1) * mrSges(3,2) + (t54 + Ifges(3,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t44 + (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t45) * t44) * qJD(1)) * qJD(1);
T = t1;
