% Calculate kinetic energy for
% S5RRRRR2
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
% Datum: 2019-12-05 18:54
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(2,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [2 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: pkin has to be [2x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:52:53
% EndTime: 2019-12-05 18:52:53
% DurationCPUTime: 0.33s
% Computational Cost: add. (236->63), mult. (415->113), div. (0->0), fcn. (250->8), ass. (0->32)
t64 = qJD(1) + qJD(2);
t67 = sin(qJ(4));
t68 = sin(qJ(3));
t71 = cos(qJ(4));
t72 = cos(qJ(3));
t58 = (t67 * t68 - t71 * t72) * t64;
t69 = sin(qJ(2));
t81 = pkin(1) * qJD(1);
t79 = t69 * t81;
t76 = qJD(3) * pkin(2) - t68 * t79;
t78 = t72 * t79;
t54 = t67 * t78 - t71 * t76;
t85 = t54 ^ 2;
t84 = t72 ^ 2;
t83 = t72 * t64;
t74 = qJD(1) ^ 2;
t82 = t74 * pkin(1) ^ 2;
t80 = t69 ^ 2 * t82;
t73 = cos(qJ(2));
t70 = cos(qJ(5));
t66 = sin(qJ(5));
t63 = qJD(3) + qJD(4);
t62 = t73 ^ 2 * t82;
t61 = -pkin(2) * t83 - t73 * t81;
t59 = (t67 * t72 + t68 * t71) * t64;
t57 = qJD(5) + t58;
t56 = t67 * t76 + t71 * t78;
t53 = t59 * t70 + t63 * t66;
t52 = -t59 * t66 + t63 * t70;
t51 = t56 * t70 + t61 * t66;
t50 = -t56 * t66 + t61 * t70;
t1 = m(6) * (t50 ^ 2 + t51 ^ 2 + t85) / 0.2e1 + t74 * Ifges(2,3) / 0.2e1 + m(3) * (t62 + t80) / 0.2e1 + m(4) * (t62 + (t68 ^ 2 + t84) * t80) / 0.2e1 + m(5) * (t56 ^ 2 + t61 ^ 2 + t85) / 0.2e1 + (-t54 * mrSges(5,1) - t56 * mrSges(5,2) + Ifges(5,3) * t63 / 0.2e1) * t63 + (t50 * mrSges(6,1) - t51 * mrSges(6,2) + Ifges(6,3) * t57 / 0.2e1) * t57 + (t61 * mrSges(5,2) + t54 * mrSges(5,3) + Ifges(5,5) * t63 + Ifges(5,1) * t59 / 0.2e1) * t59 + (t54 * mrSges(6,2) - t50 * mrSges(6,3) + Ifges(6,5) * t57 + Ifges(6,1) * t53 / 0.2e1) * t53 - (-t61 * mrSges(5,1) + t56 * mrSges(5,3) + Ifges(5,4) * t59 + Ifges(5,6) * t63 - Ifges(5,2) * t58 / 0.2e1) * t58 + (-t54 * mrSges(6,1) + t51 * mrSges(6,3) + Ifges(6,4) * t53 + Ifges(6,6) * t57 + Ifges(6,2) * t52 / 0.2e1) * t52 + (Ifges(4,2) * t84 / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(4,4) * t72 + Ifges(4,1) * t68 / 0.2e1) * t68) * t64 ^ 2 + (Ifges(4,3) * qJD(3) / 0.2e1 + (Ifges(4,5) * t68 + Ifges(4,6) * t72) * t64) * qJD(3) + ((mrSges(4,1) * t72 - mrSges(4,2) * t68 + mrSges(3,1)) * t73 * t64 + (-t64 * mrSges(3,2) - t68 * (-mrSges(4,3) * t64 * t68 + qJD(3) * mrSges(4,1)) + t72 * (-qJD(3) * mrSges(4,2) + mrSges(4,3) * t83)) * t69) * t81;
T = t1;
