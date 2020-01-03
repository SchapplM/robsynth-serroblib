% Calculate kinetic energy for
% S5RRPPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3]';
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
% Datum: 2019-12-31 19:45
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:42:58
% EndTime: 2019-12-31 19:42:58
% DurationCPUTime: 0.37s
% Computational Cost: add. (276->84), mult. (622->120), div. (0->0), fcn. (366->6), ass. (0->29)
t96 = pkin(6) * mrSges(3,3);
t86 = sin(qJ(2));
t88 = cos(qJ(2));
t73 = (-pkin(2) * t88 - qJ(3) * t86 - pkin(1)) * qJD(1);
t93 = t88 * qJD(1);
t79 = pkin(6) * t93 + qJD(2) * qJ(3);
t84 = sin(pkin(8));
t95 = cos(pkin(8));
t71 = t84 * t73 + t95 * t79;
t94 = t86 * qJD(1);
t92 = qJD(2) * pkin(2) - pkin(6) * t94 - qJD(3);
t70 = t95 * t73 - t84 * t79;
t67 = -qJ(4) * t93 + t71;
t66 = pkin(3) * t93 + qJD(4) - t70;
t75 = t84 * qJD(2) + t95 * t94;
t91 = t75 * qJ(4) + t92;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t80 = qJD(5) + t93;
t74 = -t95 * qJD(2) + t84 * t94;
t69 = t85 * t74 + t87 * t75;
t68 = t87 * t74 - t85 * t75;
t65 = t74 * pkin(3) - t91;
t64 = t74 * pkin(7) + t67;
t63 = (-pkin(3) - pkin(4)) * t74 + t91;
t62 = pkin(4) * t93 - t75 * pkin(7) + t66;
t61 = t85 * t62 + t87 * t64;
t60 = t87 * t62 - t85 * t64;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t60 ^ 2 + t61 ^ 2 + t63 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t80 / 0.2e1) * t80 + (t63 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t80 + Ifges(6,1) * t69 / 0.2e1) * t69 + (-t63 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t69 + Ifges(6,6) * t80 + Ifges(6,2) * t68 / 0.2e1) * t68 + (-t92 * mrSges(4,2) + t66 * mrSges(5,2) - t70 * mrSges(4,3) - t65 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t75) * t75 + (-t92 * mrSges(4,1) + t65 * mrSges(5,1) - t67 * mrSges(5,2) - t71 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t74 + (-Ifges(4,4) + Ifges(5,5)) * t75 + (Ifges(4,6) - Ifges(5,6)) * t93) * t74 + ((-pkin(6) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t86 + (-t70 * mrSges(4,1) + t66 * mrSges(5,1) + t71 * mrSges(4,2) - t67 * mrSges(5,3) + (-Ifges(4,5) - Ifges(5,4)) * t75 + (-pkin(6) * mrSges(3,2) + Ifges(3,6)) * qJD(2)) * t88 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t86 ^ 2 + t88 ^ 2) * pkin(6) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t96 + Ifges(3,1) / 0.2e1) * t86) * t86 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t86 + (t96 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1) * t88) * t88) * qJD(1)) * qJD(1);
T = t1;
