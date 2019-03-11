% Calculate kinetic energy for
% S6RRPRPP2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,theta3]';
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
% Datum: 2019-03-09 09:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRPP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRPP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 09:49:39
% EndTime: 2019-03-09 09:49:40
% DurationCPUTime: 0.50s
% Computational Cost: add. (482->103), mult. (1105->137), div. (0->0), fcn. (742->6), ass. (0->35)
t112 = -pkin(4) - pkin(5);
t111 = pkin(7) * mrSges(3,3);
t110 = cos(qJ(4));
t109 = pkin(7) + qJ(3);
t100 = cos(qJ(2));
t106 = t100 * qJD(1);
t99 = sin(qJ(2));
t107 = t99 * qJD(1);
t108 = cos(pkin(9));
t97 = sin(pkin(9));
t88 = t108 * t106 - t97 * t107;
t89 = (t100 * t97 + t108 * t99) * qJD(1);
t94 = qJD(3) + (-pkin(2) * t100 - pkin(1)) * qJD(1);
t75 = -pkin(3) * t88 - pkin(8) * t89 + t94;
t92 = qJD(2) * pkin(2) - t109 * t107;
t93 = t109 * t106;
t81 = t108 * t93 + t97 * t92;
t79 = qJD(2) * pkin(8) + t81;
t98 = sin(qJ(4));
t72 = t110 * t79 + t98 * t75;
t80 = t108 * t92 - t97 * t93;
t87 = qJD(4) - t88;
t70 = t87 * qJ(5) + t72;
t105 = qJD(2) * pkin(3) + t80;
t71 = t110 * t75 - t98 * t79;
t104 = qJD(5) - t71;
t83 = t98 * qJD(2) + t110 * t89;
t103 = qJ(5) * t83 + t105;
t82 = -t110 * qJD(2) + t89 * t98;
t73 = pkin(4) * t82 - t103;
t69 = -t87 * pkin(4) + t104;
t68 = t112 * t82 + qJD(6) + t103;
t67 = qJ(6) * t82 + t70;
t66 = -t83 * qJ(6) + t112 * t87 + t104;
t1 = m(4) * (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t105 ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t68 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + (t94 * mrSges(4,2) - t80 * mrSges(4,3) + Ifges(4,1) * t89 / 0.2e1) * t89 + (-t94 * mrSges(4,1) + t81 * mrSges(4,3) + Ifges(4,4) * t89 + Ifges(4,2) * t88 / 0.2e1) * t88 + (t80 * mrSges(4,1) - t81 * mrSges(4,2) + Ifges(4,5) * t89 + Ifges(4,6) * t88 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t99 + Ifges(3,6) * t100 + (-mrSges(3,1) * t99 - mrSges(3,2) * t100) * pkin(7)) * qJD(1)) * qJD(2) + (t71 * mrSges(5,1) - t69 * mrSges(6,1) - t66 * mrSges(7,1) - t72 * mrSges(5,2) + t67 * mrSges(7,2) + t70 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t87) * t87 + (-t105 * mrSges(5,2) + t69 * mrSges(6,2) + t68 * mrSges(7,2) - t71 * mrSges(5,3) - t73 * mrSges(6,3) - t66 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t83 + (Ifges(6,4) + Ifges(5,5) - Ifges(7,5)) * t87) * t83 + (-t105 * mrSges(5,1) + t73 * mrSges(6,1) - t68 * mrSges(7,1) - t70 * mrSges(6,2) - t72 * mrSges(5,3) + t67 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t82 + (-Ifges(5,6) + Ifges(6,6) - Ifges(7,6)) * t87 + (-Ifges(5,4) + Ifges(7,4) + Ifges(6,5)) * t83) * t82 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t100 ^ 2 + t99 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t111 + Ifges(3,1) / 0.2e1) * t99) * t99 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t99 + (t111 + Ifges(3,2) / 0.2e1) * t100) * t100) * qJD(1) ^ 2;
T  = t1;
