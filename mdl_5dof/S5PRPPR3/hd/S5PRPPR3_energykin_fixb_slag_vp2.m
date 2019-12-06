% Calculate kinetic energy for
% S5PRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d5,theta1,theta3]';
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
% Datum: 2019-12-05 15:27
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRPPR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:26:13
% EndTime: 2019-12-05 15:26:13
% DurationCPUTime: 0.10s
% Computational Cost: add. (87->42), mult. (177->66), div. (0->0), fcn. (80->6), ass. (0->19)
t70 = cos(qJ(2));
t63 = qJD(2) * pkin(2) + t70 * qJD(1);
t65 = sin(pkin(8));
t66 = cos(pkin(8));
t68 = sin(qJ(2));
t75 = qJD(1) * t68;
t62 = t65 * t63 + t66 * t75;
t59 = qJD(2) * qJ(4) + t62;
t76 = t59 ^ 2;
t61 = t66 * t63 - t65 * t75;
t74 = qJD(4) - t61;
t71 = qJD(3) ^ 2;
t69 = cos(qJ(5));
t67 = sin(qJ(5));
t58 = -qJD(2) * pkin(3) + t74;
t57 = (-pkin(3) - pkin(6)) * qJD(2) + t74;
t56 = t69 * qJD(3) + t67 * t57;
t55 = -t67 * qJD(3) + t69 * t57;
t1 = m(6) * (t55 ^ 2 + t56 ^ 2 + t76) / 0.2e1 + m(4) * (t61 ^ 2 + t62 ^ 2 + t71) / 0.2e1 + m(5) * (t58 ^ 2 + t71 + t76) / 0.2e1 + (m(3) * (t68 ^ 2 + t70 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t55 * mrSges(6,1) - t56 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t61 * mrSges(4,1) - t62 * mrSges(4,2) + t58 * mrSges(5,2) + t59 * mrSges(5,3) + (t70 * mrSges(3,1) - t68 * mrSges(3,2)) * qJD(1) + (t59 * mrSges(6,2) - t55 * mrSges(6,3) + Ifges(6,5) * qJD(5)) * t69 + (t59 * mrSges(6,1) - t56 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t67 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(6,1) * t69 ^ 2 / 0.2e1 + (-Ifges(6,4) * t69 + Ifges(6,2) * t67 / 0.2e1) * t67) * qJD(2)) * qJD(2);
T = t1;
