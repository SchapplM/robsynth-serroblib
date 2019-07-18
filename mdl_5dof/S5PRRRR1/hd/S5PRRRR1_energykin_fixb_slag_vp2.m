% Calculate kinetic energy for
% S5PRRRR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [2x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a4]';
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
% Datum: 2019-07-18 13:29
% Revision: 08c8d617a845f5dd194efdf9aca2774760f7818f (2019-07-16)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-07-18 13:28:17
% EndTime: 2019-07-18 13:28:18
% DurationCPUTime: 0.29s
% Computational Cost: add. (155->62), mult. (374->111), div. (0->0), fcn. (250->8), ass. (0->29)
t60 = sin(qJ(4));
t61 = sin(qJ(3));
t64 = cos(qJ(4));
t65 = cos(qJ(3));
t52 = (t60 * t61 - t64 * t65) * qJD(2);
t62 = sin(qJ(2));
t72 = qJD(1) * t62;
t68 = qJD(3) * pkin(2) - t61 * t72;
t70 = t65 * t72;
t48 = t60 * t70 - t64 * t68;
t75 = t48 ^ 2;
t74 = t65 ^ 2;
t67 = qJD(1) ^ 2;
t73 = t62 ^ 2 * t67;
t71 = qJD(2) * t65;
t66 = cos(qJ(2));
t63 = cos(qJ(5));
t59 = sin(qJ(5));
t57 = qJD(3) + qJD(4);
t56 = t66 ^ 2 * t67;
t55 = -pkin(2) * t71 - qJD(1) * t66;
t53 = (t60 * t65 + t61 * t64) * qJD(2);
t51 = qJD(5) + t52;
t50 = t60 * t68 + t64 * t70;
t47 = t53 * t63 + t57 * t59;
t46 = -t53 * t59 + t57 * t63;
t45 = t50 * t63 + t55 * t59;
t44 = -t50 * t59 + t55 * t63;
t1 = m(5) * (t50 ^ 2 + t55 ^ 2 + t75) / 0.2e1 + m(6) * (t44 ^ 2 + t45 ^ 2 + t75) / 0.2e1 + m(3) * (t56 + t73) / 0.2e1 + m(4) * (t56 + (t61 ^ 2 + t74) * t73) / 0.2e1 + m(2) * t67 / 0.2e1 + (-t48 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,3) * t57 / 0.2e1) * t57 + (t44 * mrSges(6,1) - t45 * mrSges(6,2) + Ifges(6,3) * t51 / 0.2e1) * t51 + (t55 * mrSges(5,2) + t48 * mrSges(5,3) + Ifges(5,5) * t57 + Ifges(5,1) * t53 / 0.2e1) * t53 + (t48 * mrSges(6,2) - t44 * mrSges(6,3) + Ifges(6,5) * t51 + Ifges(6,1) * t47 / 0.2e1) * t47 - (-t55 * mrSges(5,1) + t50 * mrSges(5,3) + Ifges(5,4) * t53 + Ifges(5,6) * t57 - Ifges(5,2) * t52 / 0.2e1) * t52 + (-t48 * mrSges(6,1) + t45 * mrSges(6,3) + Ifges(6,4) * t47 + Ifges(6,6) * t51 + Ifges(6,2) * t46 / 0.2e1) * t46 + (Ifges(4,2) * t74 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t65 + Ifges(4,1) * t61 / 0.2e1) * t61) * qJD(2) ^ 2 + (Ifges(4,3) * qJD(3) / 0.2e1 + (Ifges(4,5) * t61 + Ifges(4,6) * t65) * qJD(2)) * qJD(3) + ((mrSges(4,1) * t65 - mrSges(4,2) * t61 + mrSges(3,1)) * t66 * qJD(2) + (t65 * (-qJD(3) * mrSges(4,2) + mrSges(4,3) * t71) - t61 * (-mrSges(4,3) * qJD(2) * t61 + qJD(3) * mrSges(4,1)) - qJD(2) * mrSges(3,2)) * t62) * qJD(1);
T  = t1;
