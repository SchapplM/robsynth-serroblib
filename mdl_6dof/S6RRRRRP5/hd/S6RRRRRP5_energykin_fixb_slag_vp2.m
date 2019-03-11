% Calculate kinetic energy for
% S6RRRRRP5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 01:25
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 01:18:03
% EndTime: 2019-03-10 01:18:04
% DurationCPUTime: 0.59s
% Computational Cost: add. (866->107), mult. (1821->156), div. (0->0), fcn. (1324->8), ass. (0->39)
t119 = pkin(7) * mrSges(3,3);
t107 = sin(qJ(5));
t111 = cos(qJ(5));
t114 = cos(qJ(2));
t117 = t114 * qJD(1);
t105 = qJD(3) - t117;
t104 = qJD(4) + t105;
t108 = sin(qJ(4));
t112 = cos(qJ(4));
t103 = pkin(7) * t117 + qJD(2) * pkin(8);
t109 = sin(qJ(3));
t113 = cos(qJ(3));
t110 = sin(qJ(2));
t97 = (-pkin(2) * t114 - pkin(8) * t110 - pkin(1)) * qJD(1);
t92 = -t103 * t109 + t113 * t97;
t118 = t110 * qJD(1);
t99 = qJD(2) * t109 + t113 * t118;
t87 = pkin(3) * t105 - pkin(9) * t99 + t92;
t93 = t113 * t103 + t109 * t97;
t98 = qJD(2) * t113 - t109 * t118;
t89 = pkin(9) * t98 + t93;
t80 = -t108 * t89 + t112 * t87;
t91 = t108 * t98 + t112 * t99;
t76 = pkin(4) * t104 - pkin(10) * t91 + t80;
t81 = t108 * t87 + t112 * t89;
t90 = -t108 * t99 + t112 * t98;
t78 = pkin(10) * t90 + t81;
t73 = t107 * t76 + t111 * t78;
t72 = -t107 * t78 + t111 * t76;
t102 = -qJD(2) * pkin(2) + pkin(7) * t118;
t94 = -pkin(3) * t98 + t102;
t84 = -pkin(4) * t90 + t94;
t101 = qJD(5) + t104;
t83 = t107 * t90 + t111 * t91;
t82 = -t107 * t91 + t111 * t90;
t79 = -pkin(5) * t82 + qJD(6) + t84;
t71 = qJ(6) * t82 + t73;
t70 = pkin(5) * t101 - qJ(6) * t83 + t72;
t1 = m(7) * (t70 ^ 2 + t71 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t94 ^ 2) / 0.2e1 + m(4) * (t102 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + (t102 * mrSges(4,2) - t92 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t94 * mrSges(5,2) - t80 * mrSges(5,3) + Ifges(5,1) * t91 / 0.2e1) * t91 + (-t102 * mrSges(4,1) + t93 * mrSges(4,3) + Ifges(4,4) * t99 + Ifges(4,2) * t98 / 0.2e1) * t98 + (-t94 * mrSges(5,1) + t81 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,2) * t90 / 0.2e1) * t90 + (t92 * mrSges(4,1) - t93 * mrSges(4,2) + Ifges(4,5) * t99 + Ifges(4,6) * t98 + Ifges(4,3) * t105 / 0.2e1) * t105 + (t80 * mrSges(5,1) - t81 * mrSges(5,2) + Ifges(5,5) * t91 + Ifges(5,6) * t90 + Ifges(5,3) * t104 / 0.2e1) * t104 + (t84 * mrSges(6,2) + t79 * mrSges(7,2) - t72 * mrSges(6,3) - t70 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t83) * t83 + (-t84 * mrSges(6,1) - t79 * mrSges(7,1) + t73 * mrSges(6,3) + t71 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t82 + (Ifges(6,4) + Ifges(7,4)) * t83) * t82 + (t72 * mrSges(6,1) + t70 * mrSges(7,1) - t73 * mrSges(6,2) - t71 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t101 + (Ifges(6,5) + Ifges(7,5)) * t83 + (Ifges(6,6) + Ifges(7,6)) * t82) * t101 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t110 ^ 2 + t114 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t119) * t114) * t114 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t114 + (Ifges(3,1) / 0.2e1 + t119) * t110) * t110) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t114 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t110) * qJD(2)) * qJD(1);
T  = t1;
