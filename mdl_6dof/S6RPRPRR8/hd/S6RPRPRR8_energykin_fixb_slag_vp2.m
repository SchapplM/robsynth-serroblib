% Calculate kinetic energy for
% S6RPRPRR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:57:40
% EndTime: 2019-03-09 03:57:41
% DurationCPUTime: 0.40s
% Computational Cost: add. (523->94), mult. (1060->148), div. (0->0), fcn. (690->8), ass. (0->41)
t120 = qJD(1) / 0.2e1;
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t108 = sin(pkin(10));
t109 = cos(pkin(10));
t115 = cos(qJ(3));
t101 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t117 = -qJ(4) * qJD(1) + t101;
t93 = qJD(3) * pkin(3) + t115 * t117;
t112 = sin(qJ(3));
t95 = t117 * t112;
t86 = t108 * t93 + t109 * t95;
t82 = qJD(3) * pkin(8) + t86;
t118 = qJD(1) * t112;
t97 = -t108 * t115 * qJD(1) - t109 * t118;
t98 = (-t108 * t112 + t109 * t115) * qJD(1);
t107 = qJD(1) * qJ(2);
t99 = pkin(3) * t118 + qJD(4) + t107;
t87 = -pkin(4) * t97 - pkin(8) * t98 + t99;
t76 = t111 * t87 + t114 * t82;
t119 = t101 * mrSges(4,3);
t96 = qJD(5) - t97;
t85 = -t108 * t95 + t109 * t93;
t75 = -t111 * t82 + t114 * t87;
t81 = -qJD(3) * pkin(4) - t85;
t116 = qJD(1) ^ 2;
t113 = cos(qJ(6));
t110 = sin(qJ(6));
t106 = t116 * qJ(2) ^ 2;
t105 = -qJD(1) * pkin(1) + qJD(2);
t94 = qJD(6) + t96;
t89 = qJD(3) * t111 + t114 * t98;
t88 = qJD(3) * t114 - t111 * t98;
t79 = t110 * t88 + t113 * t89;
t78 = -t110 * t89 + t113 * t88;
t77 = -pkin(5) * t88 + t81;
t74 = pkin(9) * t88 + t76;
t73 = pkin(5) * t96 - pkin(9) * t89 + t75;
t72 = t110 * t73 + t113 * t74;
t71 = -t110 * t74 + t113 * t73;
t1 = m(4) * (t106 + (t112 ^ 2 + t115 ^ 2) * t101 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t99 ^ 2) / 0.2e1 + m(3) * (t105 ^ 2 + t106) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t71 ^ 2 + t72 ^ 2 + t77 ^ 2) / 0.2e1 + (t99 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t98 / 0.2e1) * t98 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t96 / 0.2e1) * t96 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,3) * t94 / 0.2e1) * t94 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t116 + (-t99 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 / 0.2e1) * t97 + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t96 + Ifges(6,1) * t89 / 0.2e1) * t89 + (t77 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,5) * t94 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t89 + Ifges(6,6) * t96 + Ifges(6,2) * t88 / 0.2e1) * t88 + (-t77 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t94 + Ifges(7,2) * t78 / 0.2e1) * t78 + (t105 * mrSges(3,2) + (mrSges(4,2) * t107 + (Ifges(4,1) * t120 - t119) * t115) * t115 + ((qJ(2) * mrSges(4,1) - Ifges(4,4) * t115) * qJD(1) + (Ifges(4,2) * t120 - t119) * t112) * t112) * qJD(1) + (t85 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (t115 * mrSges(4,1) - t112 * mrSges(4,2)) * t101 + (Ifges(4,5) * t115 - Ifges(4,6) * t112) * qJD(1)) * qJD(3);
T  = t1;
