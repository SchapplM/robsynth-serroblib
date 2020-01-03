% Calculate kinetic energy for
% S5RPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3]';
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
% Datum: 2019-12-31 18:16
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPP5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:16:05
% EndTime: 2019-12-31 18:16:05
% DurationCPUTime: 0.17s
% Computational Cost: add. (115->65), mult. (225->83), div. (0->0), fcn. (62->2), ass. (0->18)
t66 = -pkin(3) - pkin(4);
t55 = qJD(2) + (-pkin(1) - pkin(6)) * qJD(1);
t65 = t55 * mrSges(4,3);
t59 = sin(qJ(3));
t52 = qJD(3) * qJ(4) + t59 * t55;
t64 = qJ(5) * qJD(1);
t63 = qJD(1) * qJ(2);
t60 = cos(qJ(3));
t62 = qJ(4) * t60 - qJ(2);
t61 = qJD(1) ^ 2;
t58 = t61 * qJ(2) ^ 2;
t56 = -qJD(1) * pkin(1) + qJD(2);
t51 = (pkin(3) * t59 - t62) * qJD(1);
t50 = -qJD(3) * pkin(3) - t60 * t55 + qJD(4);
t49 = t59 * t64 + t52;
t48 = qJD(5) + (t66 * t59 + t62) * qJD(1);
t47 = qJD(4) + (-t55 - t64) * t60 + t66 * qJD(3);
t1 = m(4) * (t58 + (t59 ^ 2 + t60 ^ 2) * t55 ^ 2) / 0.2e1 + m(6) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + m(5) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(3) * (t56 ^ 2 + t58) / 0.2e1 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t61 + (t56 * mrSges(3,2) + (mrSges(4,2) * t63 + t50 * mrSges(5,2) + t48 * mrSges(6,2) - t51 * mrSges(5,3) - t47 * mrSges(6,3) + (-t65 + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(1)) * t60) * t60 + (t49 * mrSges(6,3) - t52 * mrSges(5,2) - t48 * mrSges(6,1) + t51 * mrSges(5,1) + mrSges(4,1) * t63 + (-t65 + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * qJD(1)) * t59 + (-Ifges(4,4) + Ifges(6,4) + Ifges(5,5)) * qJD(1) * t60) * t59) * qJD(1) + (-t50 * mrSges(5,1) - t47 * mrSges(6,1) + t49 * mrSges(6,2) + t52 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(3) + (t60 * mrSges(4,1) - t59 * mrSges(4,2)) * t55 + ((Ifges(5,4) + Ifges(4,5) - Ifges(6,5)) * t60 + (-Ifges(4,6) + Ifges(5,6) - Ifges(6,6)) * t59) * qJD(1)) * qJD(3);
T = t1;
