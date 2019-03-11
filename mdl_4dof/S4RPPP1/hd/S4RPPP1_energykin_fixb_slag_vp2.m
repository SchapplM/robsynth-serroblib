% Calculate kinetic energy for
% S4RPPP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [6x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,theta2]';
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
% Datum: 2019-03-08 18:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(6,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [6 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: pkin has to be [6x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPPP1_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 18:26:14
% EndTime: 2019-03-08 18:26:14
% DurationCPUTime: 0.20s
% Computational Cost: add. (102->64), mult. (326->82), div. (0->0), fcn. (182->4), ass. (0->24)
t60 = sin(pkin(4));
t71 = pkin(3) * t60;
t59 = sin(pkin(6));
t61 = cos(pkin(6));
t70 = t60 * qJD(1);
t64 = qJ(2) * t70;
t62 = cos(pkin(4));
t69 = t62 * qJD(1);
t67 = pkin(1) * t69;
t54 = t59 * t67 + t61 * t64;
t55 = t59 * t64;
t68 = qJD(3) + t55;
t66 = t61 * t70;
t65 = -pkin(1) * t61 - pkin(2);
t63 = -qJ(3) * t59 - pkin(1);
t58 = -pkin(1) * t70 + qJD(2);
t53 = t61 * t67 - t55;
t52 = -qJ(3) * t69 - t54;
t51 = qJD(2) + (-pkin(2) * t61 + t63) * t70;
t50 = t65 * t69 + t68;
t49 = qJD(2) + ((-pkin(2) - qJ(4)) * t61 + t63) * t70;
t48 = qJD(4) + (qJ(3) * t62 + t61 * t71) * qJD(1) + t54;
t47 = (t59 * t71 + (-qJ(4) + t65) * t62) * qJD(1) + t68;
t1 = m(3) * (t53 ^ 2 + t54 ^ 2 + t58 ^ 2) / 0.2e1 + m(4) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t47 ^ 2 + t48 ^ 2 + t49 ^ 2) / 0.2e1 + (Ifges(2,3) * qJD(1) / 0.2e1 + ((-t58 * mrSges(3,1) - t52 * mrSges(4,1) + t48 * mrSges(5,1) + t51 * mrSges(4,2) + t54 * mrSges(3,3) - t49 * mrSges(5,3) + (Ifges(3,2) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t66) * t61 + (t58 * mrSges(3,2) + t50 * mrSges(4,1) - t51 * mrSges(4,3) + t47 * mrSges(5,1) - t49 * mrSges(5,2) - t53 * mrSges(3,3) + (Ifges(3,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t59 * t70 + (Ifges(3,4) + Ifges(4,6) - Ifges(5,6)) * t66) * t59) * t60 + (t50 * mrSges(4,2) - t52 * mrSges(4,3) + t48 * mrSges(5,2) - t47 * mrSges(5,3) - t54 * mrSges(3,2) + t53 * mrSges(3,1) + (Ifges(3,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t69 + ((-Ifges(5,4) - Ifges(4,5) + Ifges(3,6)) * t61 + (-Ifges(4,4) + Ifges(3,5) + Ifges(5,5)) * t59) * t70) * t62) * qJD(1);
T  = t1;
