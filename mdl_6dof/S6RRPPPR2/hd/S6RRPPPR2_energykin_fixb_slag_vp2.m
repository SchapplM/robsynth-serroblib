% Calculate kinetic energy for
% S6RRPPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3,theta5]';
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
% Datum: 2019-03-09 08:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:09:50
% EndTime: 2019-03-09 08:09:51
% DurationCPUTime: 0.48s
% Computational Cost: add. (560->107), mult. (1305->152), div. (0->0), fcn. (892->8), ass. (0->41)
t118 = pkin(7) * mrSges(3,3);
t117 = pkin(3) + qJ(5);
t116 = pkin(7) + qJ(3);
t102 = sin(pkin(10));
t104 = cos(pkin(10));
t103 = sin(pkin(9));
t106 = sin(qJ(2));
t108 = cos(qJ(2));
t115 = cos(pkin(9));
t94 = (t103 * t108 + t106 * t115) * qJD(1);
t99 = qJD(3) + (-pkin(2) * t108 - pkin(1)) * qJD(1);
t111 = -qJ(4) * t94 + t99;
t113 = qJD(1) * t108;
t114 = qJD(1) * t106;
t93 = t103 * t114 - t113 * t115;
t79 = t117 * t93 + t111;
t97 = qJD(2) * pkin(2) - t114 * t116;
t98 = t116 * t113;
t87 = -t103 * t98 + t115 * t97;
t112 = qJD(4) - t87;
t80 = t94 * pkin(4) - qJD(2) * t117 + t112;
t74 = t102 * t80 + t104 * t79;
t88 = t103 * t97 + t115 * t98;
t86 = -qJD(2) * qJ(4) - t88;
t73 = -t102 * t79 + t104 * t80;
t83 = -pkin(4) * t93 + qJD(5) - t86;
t107 = cos(qJ(6));
t105 = sin(qJ(6));
t92 = qJD(6) + t94;
t90 = qJD(2) * t104 + t102 * t93;
t89 = -qJD(2) * t102 + t104 * t93;
t85 = -qJD(2) * pkin(3) + t112;
t84 = pkin(3) * t93 + t111;
t82 = t105 * t89 + t107 * t90;
t81 = -t105 * t90 + t107 * t89;
t75 = -pkin(5) * t89 + t83;
t72 = pkin(8) * t89 + t74;
t71 = pkin(5) * t94 - pkin(8) * t90 + t73;
t70 = t105 * t71 + t107 * t72;
t69 = -t105 * t72 + t107 * t71;
t1 = m(4) * (t87 ^ 2 + t88 ^ 2 + t99 ^ 2) / 0.2e1 + m(5) * (t84 ^ 2 + t85 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t83 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t75 ^ 2) / 0.2e1 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t92 / 0.2e1) * t92 + (t83 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,1) * t90 / 0.2e1) * t90 + (-t83 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,2) * t89 / 0.2e1) * t89 + (t75 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t92 + Ifges(7,1) * t82 / 0.2e1) * t82 + (-t75 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t82 + Ifges(7,6) * t92 + Ifges(7,2) * t81 / 0.2e1) * t81 + (t99 * mrSges(4,1) + t86 * mrSges(5,1) - t84 * mrSges(5,2) - t88 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t93) * t93 + (t87 * mrSges(4,1) - t88 * mrSges(4,2) + t85 * mrSges(5,2) - t86 * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(2) + (Ifges(5,5) - Ifges(4,6)) * t93 + (Ifges(3,5) * t106 + Ifges(3,6) * t108 + (-mrSges(3,1) * t106 - mrSges(3,2) * t108) * pkin(7)) * qJD(1)) * qJD(2) + (t85 * mrSges(5,1) + t73 * mrSges(6,1) + t99 * mrSges(4,2) - t74 * mrSges(6,2) - t87 * mrSges(4,3) - t84 * mrSges(5,3) + Ifges(6,5) * t90 + Ifges(6,6) * t89 + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t94 + (-Ifges(4,4) - Ifges(5,6)) * t93 + (-Ifges(5,4) + Ifges(4,5)) * qJD(2)) * t94 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t106 ^ 2 + t108 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t118 + Ifges(3,2) / 0.2e1) * t108) * t108 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t108 + (t118 + Ifges(3,1) / 0.2e1) * t106) * t106) * qJD(1) ^ 2;
T  = t1;
