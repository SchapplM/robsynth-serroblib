% Calculate kinetic energy for
% S6RPRRPR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6]';
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
% Datum: 2019-03-09 05:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:34:28
% EndTime: 2019-03-09 05:34:28
% DurationCPUTime: 0.35s
% Computational Cost: add. (363->92), mult. (662->133), div. (0->0), fcn. (366->6), ass. (0->36)
t109 = -pkin(4) - pkin(5);
t108 = cos(qJ(4));
t92 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t107 = t92 * mrSges(4,3);
t106 = qJD(1) / 0.2e1;
t101 = cos(qJ(3));
t99 = sin(qJ(3));
t84 = (pkin(3) * t99 - pkin(8) * t101 + qJ(2)) * qJD(1);
t85 = qJD(3) * pkin(8) + t99 * t92;
t98 = sin(qJ(4));
t78 = t108 * t85 + t98 * t84;
t105 = t101 * qJD(1);
t93 = qJD(1) * t99 + qJD(4);
t75 = t93 * qJ(5) + t78;
t77 = t108 * t84 - t98 * t85;
t86 = -qJD(3) * pkin(3) - t101 * t92;
t104 = qJD(5) - t77;
t88 = t98 * qJD(3) + t108 * t105;
t103 = qJ(5) * t88 - t86;
t102 = qJD(1) ^ 2;
t100 = cos(qJ(6));
t97 = sin(qJ(6));
t96 = t102 * qJ(2) ^ 2;
t94 = -qJD(1) * pkin(1) + qJD(2);
t91 = qJD(6) - t93;
t87 = -t108 * qJD(3) + t98 * t105;
t80 = t100 * t88 + t87 * t97;
t79 = t100 * t87 - t88 * t97;
t76 = pkin(4) * t87 - t103;
t74 = -t93 * pkin(4) + t104;
t73 = t109 * t87 + t103;
t72 = pkin(9) * t87 + t75;
t71 = -t88 * pkin(9) + t109 * t93 + t104;
t70 = t100 * t72 + t71 * t97;
t69 = t100 * t71 - t72 * t97;
t1 = m(3) * (t94 ^ 2 + t96) / 0.2e1 + m(4) * (t96 + (t101 ^ 2 + t99 ^ 2) * t92 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t91 / 0.2e1) * t91 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t102 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t101 * mrSges(4,1) - t99 * mrSges(4,2)) * t92) * qJD(3) + (t94 * mrSges(3,2) + (qJ(2) * mrSges(4,1) * qJD(1) - Ifges(4,6) * qJD(3) + (Ifges(4,2) * t106 - t107) * t99) * t99 + (Ifges(4,5) * qJD(3) + (qJ(2) * mrSges(4,2) - Ifges(4,4) * t99) * qJD(1) + (Ifges(4,1) * t106 - t107) * t101) * t101) * qJD(1) + (t73 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t91 + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t73 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,6) * t91 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t77 * mrSges(5,1) - t74 * mrSges(6,1) - t78 * mrSges(5,2) + t75 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t93) * t93 + (t86 * mrSges(5,2) + t74 * mrSges(6,2) - t77 * mrSges(5,3) - t76 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t88 + (Ifges(6,4) + Ifges(5,5)) * t93) * t88 + (t86 * mrSges(5,1) + t76 * mrSges(6,1) - t75 * mrSges(6,2) - t78 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t87 + (-Ifges(5,6) + Ifges(6,6)) * t93 + (-Ifges(5,4) + Ifges(6,5)) * t88) * t87;
T  = t1;
