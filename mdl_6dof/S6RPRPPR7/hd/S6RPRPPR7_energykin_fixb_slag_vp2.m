% Calculate kinetic energy for
% S6RPRPPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta4]';
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
% Datum: 2019-03-09 02:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:56:03
% EndTime: 2019-03-09 02:56:03
% DurationCPUTime: 0.32s
% Computational Cost: add. (343->91), mult. (694->131), div. (0->0), fcn. (396->6), ass. (0->37)
t105 = pkin(4) + pkin(8);
t86 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t104 = t86 * mrSges(4,3);
t103 = qJD(1) / 0.2e1;
t102 = sin(pkin(9));
t100 = -qJ(4) * qJD(1) + t86;
t96 = cos(qJ(3));
t79 = qJD(3) * pkin(3) + t100 * t96;
t94 = sin(qJ(3));
t80 = t100 * t94;
t92 = cos(pkin(9));
t72 = t102 * t79 + t92 * t80;
t101 = qJD(1) * t94;
t91 = qJD(1) * qJ(2);
t84 = pkin(3) * t101 + qJD(4) + t91;
t71 = -t102 * t80 + t79 * t92;
t83 = qJD(1) * t92 * t96 - t102 * t101;
t99 = qJD(5) - t71;
t70 = -qJD(3) * qJ(5) - t72;
t98 = -qJ(5) * t83 + t84;
t97 = qJD(1) ^ 2;
t95 = cos(qJ(6));
t93 = sin(qJ(6));
t90 = t97 * qJ(2) ^ 2;
t89 = -qJD(1) * pkin(1) + qJD(2);
t82 = (t102 * t96 + t92 * t94) * qJD(1);
t81 = qJD(6) + t83;
t75 = qJD(3) * t95 + t82 * t93;
t74 = -qJD(3) * t93 + t82 * t95;
t73 = pkin(4) * t82 + t98;
t69 = -qJD(3) * pkin(4) + t99;
t68 = t105 * t82 + t98;
t67 = -pkin(5) * t82 - t70;
t66 = pkin(5) * t83 - t105 * qJD(3) + t99;
t65 = t66 * t93 + t68 * t95;
t64 = t66 * t95 - t68 * t93;
t1 = m(3) * (t89 ^ 2 + t90) / 0.2e1 + m(4) * (t90 + (t94 ^ 2 + t96 ^ 2) * t86 ^ 2) / 0.2e1 + m(5) * (t71 ^ 2 + t72 ^ 2 + t84 ^ 2) / 0.2e1 + m(7) * (t64 ^ 2 + t65 ^ 2 + t67 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t97 + (t64 * mrSges(7,1) - t65 * mrSges(7,2) + Ifges(7,3) * t81 / 0.2e1) * t81 + (t67 * mrSges(7,2) - t64 * mrSges(7,3) + Ifges(7,5) * t81 + Ifges(7,1) * t75 / 0.2e1) * t75 + (-t67 * mrSges(7,1) + t65 * mrSges(7,3) + Ifges(7,4) * t75 + Ifges(7,6) * t81 + Ifges(7,2) * t74 / 0.2e1) * t74 + (t69 * mrSges(6,1) + t84 * mrSges(5,2) - t71 * mrSges(5,3) - t73 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t83) * t83 + (t84 * mrSges(5,1) + t70 * mrSges(6,1) - t73 * mrSges(6,2) - t72 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t82 + (-Ifges(5,4) - Ifges(6,6)) * t83) * t82 + (t89 * mrSges(3,2) + (mrSges(4,2) * t91 + (Ifges(4,1) * t103 - t104) * t96) * t96 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t96) * qJD(1) + (Ifges(4,2) * t103 - t104) * t94) * t94) * qJD(1) + (t71 * mrSges(5,1) - t72 * mrSges(5,2) + t69 * mrSges(6,2) - t70 * mrSges(6,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(3) + (t96 * mrSges(4,1) - t94 * mrSges(4,2)) * t86 + (-Ifges(6,4) + Ifges(5,5)) * t83 + (Ifges(6,5) - Ifges(5,6)) * t82 + (Ifges(4,5) * t96 - Ifges(4,6) * t94) * qJD(1)) * qJD(3);
T  = t1;
