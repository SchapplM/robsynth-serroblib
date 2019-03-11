% Calculate kinetic energy for
% S6RRPRPP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 10:06
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 10:03:25
% EndTime: 2019-03-09 10:03:25
% DurationCPUTime: 0.44s
% Computational Cost: add. (316->101), mult. (621->124), div. (0->0), fcn. (310->4), ass. (0->31)
t99 = -pkin(2) - pkin(8);
t98 = -pkin(4) - pkin(5);
t97 = pkin(7) * mrSges(3,3);
t88 = cos(qJ(2));
t86 = sin(qJ(2));
t93 = -qJ(3) * t86 - pkin(1);
t71 = (t99 * t88 + t93) * qJD(1);
t96 = qJD(1) * t86;
t94 = pkin(7) * t96 + qJD(3);
t72 = pkin(3) * t96 + t99 * qJD(2) + t94;
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t66 = t87 * t71 + t85 * t72;
t95 = qJD(1) * t88;
t78 = -pkin(7) * t95 - qJD(2) * qJ(3);
t80 = qJD(4) + t96;
t64 = t80 * qJ(5) + t66;
t73 = pkin(3) * t95 - t78;
t65 = -t85 * t71 + t72 * t87;
t92 = qJD(5) - t65;
t76 = qJD(2) * t87 - t85 * t95;
t91 = qJ(5) * t76 - t73;
t77 = -qJD(2) * pkin(2) + t94;
t75 = qJD(2) * t85 + t87 * t95;
t74 = (-pkin(2) * t88 + t93) * qJD(1);
t67 = pkin(4) * t75 - t91;
t63 = -pkin(4) * t80 + t92;
t62 = t98 * t75 + qJD(6) + t91;
t61 = qJ(6) * t75 + t64;
t60 = -qJ(6) * t76 + t98 * t80 + t92;
t1 = m(4) * (t74 ^ 2 + t77 ^ 2 + t78 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t73 ^ 2) / 0.2e1 + m(7) * (t60 ^ 2 + t61 ^ 2 + t62 ^ 2) / 0.2e1 + (t77 * mrSges(4,2) - t78 * mrSges(4,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(2)) * qJD(2) + (t65 * mrSges(5,1) - t63 * mrSges(6,1) - t60 * mrSges(7,1) - t66 * mrSges(5,2) + t61 * mrSges(7,2) + t64 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t80) * t80 + (t73 * mrSges(5,2) + t63 * mrSges(6,2) + t62 * mrSges(7,2) - t65 * mrSges(5,3) - t67 * mrSges(6,3) - t60 * mrSges(7,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t76 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t80) * t76 + (t73 * mrSges(5,1) + t67 * mrSges(6,1) - t62 * mrSges(7,1) - t64 * mrSges(6,2) - t66 * mrSges(5,3) + t61 * mrSges(7,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t75 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t80 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t76) * t75 + ((-t78 * mrSges(4,1) + t74 * mrSges(4,2) + (-mrSges(3,2) * pkin(7) - Ifges(4,5) + Ifges(3,6)) * qJD(2)) * t88 + (t77 * mrSges(4,1) - t74 * mrSges(4,3) + (-mrSges(3,1) * pkin(7) - Ifges(4,4) + Ifges(3,5)) * qJD(2)) * t86 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t86 ^ 2 + t88 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(4,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + t97) * t88) * t88 + (-pkin(1) * mrSges(3,2) + (Ifges(4,2) / 0.2e1 + Ifges(3,1) / 0.2e1 + t97) * t86 + (Ifges(3,4) + Ifges(4,6)) * t88) * t86) * qJD(1)) * qJD(1);
T  = t1;
