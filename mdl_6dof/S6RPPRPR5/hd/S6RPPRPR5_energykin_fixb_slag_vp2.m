% Calculate kinetic energy for
% S6RPPRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta5]';
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
% Datum: 2019-03-09 01:49
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:48:27
% EndTime: 2019-03-09 01:48:27
% DurationCPUTime: 0.25s
% Computational Cost: add. (278->79), mult. (513->121), div. (0->0), fcn. (266->6), ass. (0->32)
t99 = -pkin(1) - qJ(3);
t84 = -qJD(1) * t99 - qJD(2);
t102 = t84 ^ 2;
t101 = m(3) / 0.2e1;
t87 = qJD(1) * qJ(2) + qJD(3);
t83 = -qJD(1) * pkin(7) + t87;
t100 = t83 * mrSges(5,3);
t93 = sin(qJ(4));
t95 = cos(qJ(4));
t76 = -qJD(2) + (pkin(4) * t93 - qJ(5) * t95 - t99) * qJD(1);
t79 = qJD(4) * qJ(5) + t93 * t83;
t90 = sin(pkin(9));
t91 = cos(pkin(9));
t70 = t90 * t76 + t91 * t79;
t98 = qJD(1) * t95;
t97 = t93 * qJD(1);
t69 = t91 * t76 - t90 * t79;
t78 = -qJD(4) * pkin(4) - t95 * t83 + qJD(5);
t94 = cos(qJ(6));
t92 = sin(qJ(6));
t88 = -qJD(1) * pkin(1) + qJD(2);
t86 = qJD(6) + t97;
t81 = t90 * qJD(4) + t91 * t98;
t80 = t91 * qJD(4) - t90 * t98;
t73 = -t80 * pkin(5) + t78;
t72 = t92 * t80 + t94 * t81;
t71 = t94 * t80 - t92 * t81;
t68 = t80 * pkin(8) + t70;
t67 = pkin(5) * t97 - t81 * pkin(8) + t69;
t66 = t92 * t67 + t94 * t68;
t65 = t94 * t67 - t92 * t68;
t1 = m(5) * (t102 + (t93 ^ 2 + t95 ^ 2) * t83 ^ 2) / 0.2e1 + t88 ^ 2 * t101 + m(4) * (t87 ^ 2 + t102) / 0.2e1 + m(7) * (t65 ^ 2 + t66 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t78 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + (qJ(2) * t101 + mrSges(3,3)) * qJ(2)) * qJD(1) ^ 2 + (t65 * mrSges(7,1) - t66 * mrSges(7,2) + Ifges(7,3) * t86 / 0.2e1) * t86 + (t78 * mrSges(6,2) - t69 * mrSges(6,3) + Ifges(6,1) * t81 / 0.2e1) * t81 + (Ifges(5,3) * qJD(4) / 0.2e1 + (t95 * mrSges(5,1) - t93 * mrSges(5,2)) * t83) * qJD(4) + (t73 * mrSges(7,2) - t65 * mrSges(7,3) + Ifges(7,5) * t86 + Ifges(7,1) * t72 / 0.2e1) * t72 + (t88 * mrSges(3,2) + t87 * mrSges(4,2) + t84 * mrSges(4,3) + (t84 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (-t100 + Ifges(5,1) * qJD(1) / 0.2e1) * t95) * t95 + (-Ifges(5,4) * t98 + t84 * mrSges(5,1) + t69 * mrSges(6,1) - t70 * mrSges(6,2) + Ifges(6,5) * t81 - Ifges(5,6) * qJD(4) + (-t100 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(1)) * t93) * t93) * qJD(1) + (Ifges(6,6) * t97 - t78 * mrSges(6,1) + t70 * mrSges(6,3) + Ifges(6,4) * t81 + Ifges(6,2) * t80 / 0.2e1) * t80 + (-t73 * mrSges(7,1) + t66 * mrSges(7,3) + Ifges(7,4) * t72 + Ifges(7,6) * t86 + Ifges(7,2) * t71 / 0.2e1) * t71;
T  = t1;
