% Calculate kinetic energy for
% S6RPRRRR9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,d6]';
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
% Datum: 2019-03-09 07:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:23:02
% EndTime: 2019-03-09 07:23:03
% DurationCPUTime: 0.48s
% Computational Cost: add. (589->95), mult. (1128->150), div. (0->0), fcn. (740->8), ass. (0->40)
t118 = qJD(1) / 0.2e1;
t108 = sin(qJ(5));
t112 = cos(qJ(5));
t110 = sin(qJ(3));
t103 = t110 * qJD(1) + qJD(4);
t109 = sin(qJ(4));
t113 = cos(qJ(4));
t114 = cos(qJ(3));
t94 = (pkin(3) * t110 - pkin(8) * t114 + qJ(2)) * qJD(1);
t102 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t95 = qJD(3) * pkin(8) + t110 * t102;
t86 = -t109 * t95 + t113 * t94;
t116 = t114 * qJD(1);
t98 = qJD(3) * t109 + t113 * t116;
t83 = pkin(4) * t103 - pkin(9) * t98 + t86;
t87 = t109 * t94 + t113 * t95;
t97 = qJD(3) * t113 - t109 * t116;
t85 = pkin(9) * t97 + t87;
t77 = t108 * t83 + t112 * t85;
t117 = t102 * mrSges(4,3);
t101 = qJD(5) + t103;
t76 = -t108 * t85 + t112 * t83;
t96 = -qJD(3) * pkin(3) - t114 * t102;
t90 = -pkin(4) * t97 + t96;
t115 = qJD(1) ^ 2;
t111 = cos(qJ(6));
t107 = sin(qJ(6));
t106 = t115 * qJ(2) ^ 2;
t104 = -qJD(1) * pkin(1) + qJD(2);
t99 = qJD(6) + t101;
t89 = t108 * t97 + t112 * t98;
t88 = -t108 * t98 + t112 * t97;
t80 = -pkin(5) * t88 + t90;
t79 = t107 * t88 + t111 * t89;
t78 = -t107 * t89 + t111 * t88;
t75 = pkin(10) * t88 + t77;
t74 = pkin(5) * t101 - pkin(10) * t89 + t76;
t73 = t107 * t74 + t111 * t75;
t72 = -t107 * t75 + t111 * t74;
t1 = m(3) * (t104 ^ 2 + t106) / 0.2e1 + m(4) * (t106 + (t110 ^ 2 + t114 ^ 2) * t102 ^ 2) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t72 ^ 2 + t73 ^ 2 + t80 ^ 2) / 0.2e1 + (t72 * mrSges(7,1) - t73 * mrSges(7,2) + Ifges(7,3) * t99 / 0.2e1) * t99 + (t96 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,1) * t98 / 0.2e1) * t98 + (t90 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,1) * t89 / 0.2e1) * t89 + (qJ(2) * mrSges(3,3) + Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1) * t115 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t114 * mrSges(4,1) - t110 * mrSges(4,2)) * t102) * qJD(3) + (t104 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t118 - t117) * t114) * t114 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t114) * qJD(1) + (Ifges(4,2) * t118 - t117) * t110) * t110) * qJD(1) + (-t96 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,4) * t98 + Ifges(5,2) * t97 / 0.2e1) * t97 + (-t90 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t89 + Ifges(6,2) * t88 / 0.2e1) * t88 + (t80 * mrSges(7,2) - t72 * mrSges(7,3) + Ifges(7,5) * t99 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t80 * mrSges(7,1) + t73 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t99 + Ifges(7,2) * t78 / 0.2e1) * t78 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,5) * t98 + Ifges(5,6) * t97 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,5) * t89 + Ifges(6,6) * t88 + Ifges(6,3) * t101 / 0.2e1) * t101;
T  = t1;
