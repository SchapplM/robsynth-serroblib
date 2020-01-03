% Calculate kinetic energy for
% S5RRRPP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3]';
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
% Datum: 2019-12-31 21:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(7,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPP8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:07:26
% EndTime: 2019-12-31 21:07:27
% DurationCPUTime: 0.30s
% Computational Cost: add. (234->80), mult. (486->104), div. (0->0), fcn. (260->4), ass. (0->25)
t80 = pkin(6) * mrSges(3,3);
t79 = cos(qJ(3));
t78 = pkin(3) + qJ(5);
t70 = sin(qJ(2));
t71 = cos(qJ(2));
t59 = (-pkin(2) * t71 - pkin(7) * t70 - pkin(1)) * qJD(1);
t76 = t71 * qJD(1);
t65 = pkin(6) * t76 + qJD(2) * pkin(7);
t69 = sin(qJ(3));
t57 = t69 * t59 + t79 * t65;
t77 = t70 * qJD(1);
t64 = -qJD(2) * pkin(2) + pkin(6) * t77;
t66 = qJD(3) - t76;
t54 = -t66 * qJ(4) - t57;
t56 = t79 * t59 - t69 * t65;
t75 = qJD(4) - t56;
t61 = t69 * qJD(2) + t79 * t77;
t74 = -t61 * qJ(4) + t64;
t60 = -t79 * qJD(2) + t69 * t77;
t55 = t60 * pkin(3) + t74;
t53 = -t66 * pkin(3) + t75;
t52 = t78 * t60 + t74;
t51 = -t60 * pkin(4) + qJD(5) - t54;
t50 = t61 * pkin(4) - t78 * t66 + t75;
t1 = m(6) * (t50 ^ 2 + t51 ^ 2 + t52 ^ 2) / 0.2e1 + m(5) * (t53 ^ 2 + t54 ^ 2 + t55 ^ 2) / 0.2e1 + m(4) * (t56 ^ 2 + t57 ^ 2 + t64 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t56 * mrSges(4,1) - t57 * mrSges(4,2) + t53 * mrSges(5,2) + t51 * mrSges(6,2) - t54 * mrSges(5,3) - t50 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t66) * t66 + (t53 * mrSges(5,1) + t50 * mrSges(6,1) + t64 * mrSges(4,2) - t52 * mrSges(6,2) - t56 * mrSges(4,3) - t55 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t61 + (-Ifges(5,4) + Ifges(4,5) + Ifges(6,5)) * t66) * t61 + (t64 * mrSges(4,1) + t54 * mrSges(5,1) - t51 * mrSges(6,1) - t55 * mrSges(5,2) - t57 * mrSges(4,3) + t52 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t60 + (Ifges(6,4) + Ifges(5,5) - Ifges(4,6)) * t66 + (-Ifges(4,4) - Ifges(5,6) + Ifges(6,6)) * t61) * t60 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t70 ^ 2 + t71 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t80 + Ifges(3,2) / 0.2e1) * t71) * t71 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t71 + (t80 + Ifges(3,1) / 0.2e1) * t70) * t70) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t71 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t70) * qJD(2)) * qJD(1);
T = t1;
