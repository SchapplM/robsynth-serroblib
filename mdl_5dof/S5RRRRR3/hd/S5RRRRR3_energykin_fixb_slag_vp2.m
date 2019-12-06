% Calculate kinetic energy for
% S5RRRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [5x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a3,a4,a5,d1,d4]';
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
% Datum: 2019-12-05 18:57
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(5,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [5 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: pkin has to be [5x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:55:02
% EndTime: 2019-12-05 18:55:03
% DurationCPUTime: 0.32s
% Computational Cost: add. (347->70), mult. (725->130), div. (0->0), fcn. (534->8), ass. (0->30)
t88 = pkin(1) * qJD(2);
t83 = cos(qJ(2));
t87 = t83 * qJD(1);
t78 = sin(qJ(3));
t79 = sin(qJ(2));
t82 = cos(qJ(3));
t70 = -t78 * t79 * qJD(1) + t82 * t87;
t71 = (t78 * t83 + t79 * t82) * qJD(1);
t64 = -pkin(1) * t87 - t70 * pkin(2) - t71 * pkin(5);
t75 = qJD(2) + qJD(3);
t72 = t75 * pkin(5) + t78 * t88;
t77 = sin(qJ(4));
t81 = cos(qJ(4));
t61 = t81 * t64 - t77 * t72;
t73 = -t75 * pkin(2) - t82 * t88;
t69 = qJD(4) - t70;
t85 = qJD(1) ^ 2;
t80 = cos(qJ(5));
t76 = sin(qJ(5));
t68 = qJD(5) + t69;
t67 = t81 * t71 + t77 * t75;
t66 = -t77 * t71 + t81 * t75;
t65 = -t66 * pkin(3) + t73;
t62 = t77 * t64 + t81 * t72;
t60 = t76 * t66 + t80 * t67;
t59 = t80 * t66 - t76 * t67;
t58 = t69 * pkin(3) + t61;
t57 = t76 * t58 + t80 * t62;
t56 = t80 * t58 - t76 * t62;
t1 = Ifges(4,3) * t75 ^ 2 / 0.2e1 + t85 * Ifges(2,3) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t65 ^ 2) / 0.2e1 + (Ifges(4,5) * t75 + Ifges(4,1) * t71 / 0.2e1) * t71 + m(4) * (t83 ^ 2 * t85 + (t78 ^ 2 + t82 ^ 2) * qJD(2) ^ 2) * pkin(1) ^ 2 / 0.2e1 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,3) * t69 / 0.2e1) * t69 + (t56 * mrSges(6,1) - t57 * mrSges(6,2) + Ifges(6,3) * t68 / 0.2e1) * t68 + (t73 * mrSges(5,2) - t61 * mrSges(5,3) + Ifges(5,5) * t69 + Ifges(5,1) * t67 / 0.2e1) * t67 + (t65 * mrSges(6,2) - t56 * mrSges(6,3) + Ifges(6,5) * t68 + Ifges(6,1) * t60 / 0.2e1) * t60 + (Ifges(4,4) * t71 + Ifges(4,6) * t75 + Ifges(4,2) * t70 / 0.2e1) * t70 + (-t73 * mrSges(5,1) + t62 * mrSges(5,3) + Ifges(5,4) * t67 + Ifges(5,6) * t69 + Ifges(5,2) * t66 / 0.2e1) * t66 + (-t65 * mrSges(6,1) + t57 * mrSges(6,3) + Ifges(6,4) * t60 + Ifges(6,6) * t68 + Ifges(6,2) * t59 / 0.2e1) * t59 + (Ifges(3,1) * t79 ^ 2 * qJD(1) / 0.2e1 + (-pkin(1) * (-t70 * mrSges(4,1) + t71 * mrSges(4,2)) + (Ifges(3,4) * t79 + Ifges(3,2) * t83 / 0.2e1) * qJD(1)) * t83) * qJD(1) + (Ifges(3,3) * qJD(2) / 0.2e1 + (Ifges(3,5) * t79 + Ifges(3,6) * t83) * qJD(1) + (t78 * (-t75 * mrSges(4,2) + t70 * mrSges(4,3)) + t82 * (t75 * mrSges(4,1) - t71 * mrSges(4,3))) * pkin(1)) * qJD(2);
T = t1;
