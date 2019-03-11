% Calculate kinetic energy for
% S6PRPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1]';
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
% Datum: 2019-03-08 19:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:51:06
% EndTime: 2019-03-08 19:51:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (200->83), mult. (396->116), div. (0->0), fcn. (200->8), ass. (0->37)
t101 = sin(qJ(2));
t97 = sin(pkin(6));
t114 = qJD(1) * t97;
t110 = t101 * t114;
t89 = qJD(2) * qJ(3) + t110;
t115 = t89 ^ 2;
t100 = sin(qJ(4));
t103 = cos(qJ(4));
t98 = cos(pkin(6));
t113 = qJD(1) * t98;
t104 = cos(qJ(2));
t107 = -t104 * t114 + qJD(3);
t85 = (-pkin(2) - pkin(8)) * qJD(2) + t107;
t82 = t100 * t85 + t103 * t113;
t112 = qJD(2) * t100;
t111 = t103 * qJD(2);
t91 = t100 * t113;
t81 = t103 * t85 - t91;
t109 = -qJ(5) * t103 + qJ(3);
t108 = pkin(4) * t112 + t110;
t79 = -qJD(4) * qJ(5) - t82;
t106 = qJD(1) ^ 2;
t102 = cos(qJ(6));
t99 = sin(qJ(6));
t95 = t98 ^ 2 * t106;
t93 = qJD(6) + t111;
t88 = t102 * qJD(4) + t112 * t99;
t87 = -t99 * qJD(4) + t102 * t112;
t86 = -qJD(2) * pkin(2) + t107;
t83 = qJD(2) * t109 + t108;
t80 = (pkin(9) * t100 + t109) * qJD(2) + t108;
t78 = -qJD(4) * pkin(4) + qJD(5) - t81;
t77 = -pkin(5) * t112 - t79;
t76 = qJD(5) + t91 + (pkin(5) * qJD(2) - t85) * t103 + (-pkin(4) - pkin(9)) * qJD(4);
t75 = t102 * t80 + t99 * t76;
t74 = t102 * t76 - t99 * t80;
t1 = m(3) * (t95 + (t101 ^ 2 + t104 ^ 2) * t97 ^ 2 * t106) / 0.2e1 + m(2) * t106 / 0.2e1 + m(4) * (t86 ^ 2 + t115 + t95) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t83 ^ 2) / 0.2e1 + m(5) * (t81 ^ 2 + t82 ^ 2 + t115) / 0.2e1 + m(7) * (t74 ^ 2 + t75 ^ 2 + t77 ^ 2) / 0.2e1 + (t74 * mrSges(7,1) - t75 * mrSges(7,2) + Ifges(7,3) * t93 / 0.2e1) * t93 + (t77 * mrSges(7,2) - t74 * mrSges(7,3) + Ifges(7,5) * t93 + Ifges(7,1) * t88 / 0.2e1) * t88 + (-t77 * mrSges(7,1) + t75 * mrSges(7,3) + Ifges(7,4) * t88 + Ifges(7,6) * t93 + Ifges(7,2) * t87 / 0.2e1) * t87 + (t81 * mrSges(5,1) - t82 * mrSges(5,2) + t78 * mrSges(6,2) - t79 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4)) * qJD(4) + (t86 * mrSges(4,2) + t89 * mrSges(4,3) + (mrSges(3,1) * t104 - mrSges(3,2) * t101) * t114 + (t78 * mrSges(6,1) + t89 * mrSges(5,2) - t81 * mrSges(5,3) - t83 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t111 + (-Ifges(6,4) + Ifges(5,5)) * qJD(4)) * t103 + (t89 * mrSges(5,1) + t79 * mrSges(6,1) - t83 * mrSges(6,2) - t82 * mrSges(5,3) + (Ifges(6,5) - Ifges(5,6)) * qJD(4)) * t100 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + ((Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t100 + (-Ifges(5,4) - Ifges(6,6)) * t103) * t100) * qJD(2)) * qJD(2);
T  = t1;
