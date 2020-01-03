% Calculate kinetic energy for
% S4RPRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [7x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 16:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RPRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(7,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [7 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: pkin has to be [7x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RPRR7_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 16:53:40
% EndTime: 2019-12-31 16:53:40
% DurationCPUTime: 0.25s
% Computational Cost: add. (165->56), mult. (421->94), div. (0->0), fcn. (262->6), ass. (0->26)
t70 = cos(pkin(7));
t82 = t70 ^ 2;
t81 = qJD(1) * (pkin(5) + qJ(2));
t69 = sin(pkin(7));
t72 = sin(qJ(3));
t74 = cos(qJ(3));
t61 = (t69 * t72 - t70 * t74) * qJD(1);
t80 = m(3) / 0.2e1;
t63 = t69 * t81;
t64 = t70 * t81;
t55 = -t72 * t63 + t74 * t64;
t54 = -t63 * t74 - t64 * t72;
t65 = qJD(2) + (-pkin(2) * t70 - pkin(1)) * qJD(1);
t73 = cos(qJ(4));
t71 = sin(qJ(4));
t66 = -qJD(1) * pkin(1) + qJD(2);
t62 = (t69 * t74 + t70 * t72) * qJD(1);
t58 = qJD(4) + t61;
t57 = qJD(3) * t71 + t62 * t73;
t56 = qJD(3) * t73 - t62 * t71;
t53 = qJD(3) * pkin(6) + t55;
t52 = -qJD(3) * pkin(3) - t54;
t51 = pkin(3) * t61 - pkin(6) * t62 + t65;
t50 = t51 * t71 + t53 * t73;
t49 = t51 * t73 - t53 * t71;
t1 = m(5) * (t49 ^ 2 + t50 ^ 2 + t52 ^ 2) / 0.2e1 + m(4) * (t54 ^ 2 + t55 ^ 2 + t65 ^ 2) / 0.2e1 + t66 ^ 2 * t80 + (t65 * mrSges(4,2) - t54 * mrSges(4,3) + Ifges(4,1) * t62 / 0.2e1) * t62 + (t49 * mrSges(5,1) - t50 * mrSges(5,2) + Ifges(5,3) * t58 / 0.2e1) * t58 - (-t65 * mrSges(4,1) + t55 * mrSges(4,3) + Ifges(4,4) * t62 - Ifges(4,2) * t61 / 0.2e1) * t61 + (t52 * mrSges(5,2) - t49 * mrSges(5,3) + Ifges(5,5) * t58 + Ifges(5,1) * t57 / 0.2e1) * t57 + (-t52 * mrSges(5,1) + t50 * mrSges(5,3) + Ifges(5,4) * t57 + Ifges(5,6) * t58 + Ifges(5,2) * t56 / 0.2e1) * t56 + (t54 * mrSges(4,1) - t55 * mrSges(4,2) + Ifges(4,5) * t62 - Ifges(4,6) * t61 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t66 * (-mrSges(3,1) * t70 + mrSges(3,2) * t69) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t80 + mrSges(3,3)) * (t69 ^ 2 + t82) * qJ(2) + Ifges(3,2) * t82 / 0.2e1 + (Ifges(3,4) * t70 + Ifges(3,1) * t69 / 0.2e1) * t69) * qJD(1)) * qJD(1);
T = t1;
