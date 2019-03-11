% Calculate kinetic energy for
% S6RRRRPP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4]';
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
% Datum: 2019-03-09 20:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRPP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRPP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 20:53:49
% EndTime: 2019-03-09 20:53:49
% DurationCPUTime: 0.48s
% Computational Cost: add. (542->103), mult. (1105->138), div. (0->0), fcn. (742->6), ass. (0->34)
t107 = qJD(1) * (-pkin(8) - pkin(7));
t94 = sin(qJ(3));
t95 = sin(qJ(2));
t96 = cos(qJ(3));
t97 = cos(qJ(2));
t84 = (t94 * t95 - t96 * t97) * qJD(1);
t105 = pkin(7) * mrSges(3,3);
t104 = cos(qJ(4));
t103 = pkin(4) + qJ(6);
t85 = (t94 * t97 + t95 * t96) * qJD(1);
t90 = (-pkin(2) * t97 - pkin(1)) * qJD(1);
t73 = pkin(3) * t84 - pkin(9) * t85 + t90;
t88 = qJD(2) * pkin(2) + t95 * t107;
t89 = t97 * t107;
t79 = t94 * t88 - t96 * t89;
t92 = qJD(2) + qJD(3);
t77 = pkin(9) * t92 + t79;
t93 = sin(qJ(4));
t71 = t104 * t77 + t93 * t73;
t78 = t88 * t96 + t94 * t89;
t83 = qJD(4) + t84;
t68 = -qJ(5) * t83 - t71;
t70 = t104 * t73 - t93 * t77;
t76 = -pkin(3) * t92 - t78;
t101 = qJD(5) - t70;
t81 = t104 * t85 + t93 * t92;
t100 = -qJ(5) * t81 + t76;
t80 = -t104 * t92 + t85 * t93;
t69 = pkin(4) * t80 + t100;
t67 = -t83 * pkin(4) + t101;
t66 = t103 * t80 + t100;
t65 = -pkin(5) * t80 + qJD(6) - t68;
t64 = t81 * pkin(5) - t103 * t83 + t101;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t78 ^ 2 + t79 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t76 ^ 2) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t66 ^ 2) / 0.2e1 + m(6) * (t67 ^ 2 + t68 ^ 2 + t69 ^ 2) / 0.2e1 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,3) * t92 / 0.2e1) * t92 + (t90 * mrSges(4,2) - t78 * mrSges(4,3) + Ifges(4,5) * t92 + Ifges(4,1) * t85 / 0.2e1) * t85 - (-t90 * mrSges(4,1) + t79 * mrSges(4,3) + Ifges(4,4) * t85 + Ifges(4,6) * t92 - Ifges(4,2) * t84 / 0.2e1) * t84 + (t70 * mrSges(5,1) - t71 * mrSges(5,2) + t67 * mrSges(6,2) + t65 * mrSges(7,2) - t68 * mrSges(6,3) - t64 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t83) * t83 + (t67 * mrSges(6,1) + t64 * mrSges(7,1) + t76 * mrSges(5,2) - t66 * mrSges(7,2) - t70 * mrSges(5,3) - t69 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t81 + (-Ifges(6,4) + Ifges(5,5) + Ifges(7,5)) * t83) * t81 + (t76 * mrSges(5,1) + t68 * mrSges(6,1) - t65 * mrSges(7,1) - t69 * mrSges(6,2) - t71 * mrSges(5,3) + t66 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t80 + (Ifges(7,4) + Ifges(6,5) - Ifges(5,6)) * t83 + (-Ifges(5,4) - Ifges(6,6) + Ifges(7,6)) * t81) * t80 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t95 ^ 2 + t97 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t105) * t97) * t97 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t97 + (Ifges(3,1) / 0.2e1 + t105) * t95) * t95) * qJD(1) + ((-mrSges(3,2) * pkin(7) + Ifges(3,6)) * t97 + (-mrSges(3,1) * pkin(7) + Ifges(3,5)) * t95) * qJD(2)) * qJD(1);
T  = t1;
