% Calculate kinetic energy for
% S6RPRRRR8
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
% Datum: 2019-03-09 07:22
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 07:19:36
% EndTime: 2019-03-09 07:19:36
% DurationCPUTime: 0.47s
% Computational Cost: add. (577->94), mult. (1060->149), div. (0->0), fcn. (690->8), ass. (0->42)
t122 = qJD(1) / 0.2e1;
t111 = sin(qJ(5));
t115 = cos(qJ(5));
t107 = qJD(3) + qJD(4);
t112 = sin(qJ(4));
t116 = cos(qJ(4));
t117 = cos(qJ(3));
t102 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t119 = -pkin(8) * qJD(1) + t102;
t94 = qJD(3) * pkin(3) + t117 * t119;
t113 = sin(qJ(3));
t96 = t119 * t113;
t87 = t112 * t94 + t116 * t96;
t83 = pkin(9) * t107 + t87;
t109 = qJD(1) * qJ(2);
t120 = qJD(1) * t113;
t100 = pkin(3) * t120 + t109;
t98 = -t112 * t117 * qJD(1) - t116 * t120;
t99 = (-t112 * t113 + t116 * t117) * qJD(1);
t88 = -pkin(4) * t98 - pkin(9) * t99 + t100;
t77 = t111 * t88 + t115 * t83;
t121 = t102 * mrSges(4,3);
t97 = qJD(5) - t98;
t76 = -t111 * t83 + t115 * t88;
t86 = -t112 * t96 + t116 * t94;
t82 = -pkin(4) * t107 - t86;
t118 = qJD(1) ^ 2;
t114 = cos(qJ(6));
t110 = sin(qJ(6));
t108 = t118 * qJ(2) ^ 2;
t106 = -qJD(1) * pkin(1) + qJD(2);
t95 = qJD(6) + t97;
t90 = t107 * t111 + t115 * t99;
t89 = t107 * t115 - t111 * t99;
t80 = t110 * t89 + t114 * t90;
t79 = -t110 * t90 + t114 * t89;
t78 = -pkin(5) * t89 + t82;
t75 = pkin(10) * t89 + t77;
t74 = pkin(5) * t97 - pkin(10) * t90 + t76;
t73 = t110 * t74 + t114 * t75;
t72 = -t110 * t75 + t114 * t74;
t1 = m(4) * (t108 + (t113 ^ 2 + t117 ^ 2) * t102 ^ 2) / 0.2e1 + m(5) * (t100 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(3) * (t106 ^ 2 + t108) / 0.2e1 + m(6) * (t76 ^ 2 + t77 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t72 ^ 2 + t73 ^ 2 + t78 ^ 2) / 0.2e1 + (t100 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,1) * t99 / 0.2e1) * t99 + (t76 * mrSges(6,1) - t77 * mrSges(6,2) + Ifges(6,3) * t97 / 0.2e1) * t97 + (t72 * mrSges(7,1) - t73 * mrSges(7,2) + Ifges(7,3) * t95 / 0.2e1) * t95 + (Ifges(3,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + qJ(2) * mrSges(3,3)) * t118 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t117 * mrSges(4,1) - t113 * mrSges(4,2)) * t102) * qJD(3) + (t106 * mrSges(3,2) + (mrSges(4,2) * t109 + Ifges(4,5) * qJD(3) + (Ifges(4,1) * t122 - t121) * t117) * t117 + (-Ifges(4,6) * qJD(3) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t117) * qJD(1) + (Ifges(4,2) * t122 - t121) * t113) * t113) * qJD(1) + (-t100 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,4) * t99 + Ifges(5,2) * t98 / 0.2e1) * t98 + (t82 * mrSges(6,2) - t76 * mrSges(6,3) + Ifges(6,5) * t97 + Ifges(6,1) * t90 / 0.2e1) * t90 + (t78 * mrSges(7,2) - t72 * mrSges(7,3) + Ifges(7,5) * t95 + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t82 * mrSges(6,1) + t77 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,6) * t97 + Ifges(6,2) * t89 / 0.2e1) * t89 + (-t78 * mrSges(7,1) + t73 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,6) * t95 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,5) * t99 + Ifges(5,6) * t98 + Ifges(5,3) * t107 / 0.2e1) * t107;
T  = t1;
