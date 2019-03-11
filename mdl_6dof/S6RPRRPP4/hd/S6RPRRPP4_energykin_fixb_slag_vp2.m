% Calculate kinetic energy for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
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
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:38:36
% EndTime: 2019-03-09 04:38:36
% DurationCPUTime: 0.55s
% Computational Cost: add. (671->100), mult. (1599->145), div. (0->0), fcn. (1182->8), ass. (0->38)
t107 = cos(pkin(9));
t120 = t107 ^ 2;
t119 = qJD(1) * (pkin(7) + qJ(2));
t106 = sin(pkin(9));
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t97 = (t106 * t109 - t107 * t111) * qJD(1);
t118 = m(3) / 0.2e1;
t105 = sin(pkin(10));
t116 = cos(pkin(10));
t108 = sin(qJ(4));
t110 = cos(qJ(4));
t101 = qJD(2) + (-pkin(2) * t107 - pkin(1)) * qJD(1);
t98 = (t106 * t111 + t107 * t109) * qJD(1);
t84 = pkin(3) * t97 - pkin(8) * t98 + t101;
t100 = t107 * t119;
t99 = t106 * t119;
t89 = t111 * t100 - t109 * t99;
t87 = qJD(3) * pkin(8) + t89;
t77 = -t108 * t87 + t110 * t84;
t92 = qJD(3) * t108 + t110 * t98;
t93 = qJD(4) + t97;
t74 = pkin(4) * t93 - qJ(5) * t92 + t77;
t78 = t108 * t84 + t110 * t87;
t91 = qJD(3) * t110 - t108 * t98;
t76 = qJ(5) * t91 + t78;
t71 = t105 * t74 + t116 * t76;
t88 = -t109 * t100 - t111 * t99;
t86 = -qJD(3) * pkin(3) - t88;
t70 = -t105 * t76 + t116 * t74;
t79 = -pkin(4) * t91 + qJD(5) + t86;
t102 = -qJD(1) * pkin(1) + qJD(2);
t81 = t105 * t91 + t116 * t92;
t80 = t105 * t92 - t116 * t91;
t72 = pkin(5) * t80 - qJ(6) * t81 + t79;
t69 = qJ(6) * t93 + t71;
t68 = -t93 * pkin(5) + qJD(6) - t70;
t1 = m(4) * (t101 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t102 ^ 2 * t118 + m(5) * (t77 ^ 2 + t78 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t79 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t88 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t86 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 - (-t101 * mrSges(4,1) + t89 * mrSges(4,3) + Ifges(4,4) * t98 - Ifges(4,2) * t97 / 0.2e1) * t97 + (-t86 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t92 + Ifges(5,2) * t91 / 0.2e1) * t91 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,5) * t98 - Ifges(4,6) * t97 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t79 * mrSges(6,2) + t68 * mrSges(7,2) - t70 * mrSges(6,3) - t72 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81) * t81 + (t79 * mrSges(6,1) + t72 * mrSges(7,1) - t69 * mrSges(7,2) - t71 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + (t77 * mrSges(5,1) + t70 * mrSges(6,1) - t68 * mrSges(7,1) - t78 * mrSges(5,2) - t71 * mrSges(6,2) + t69 * mrSges(7,3) + Ifges(5,5) * t92 + Ifges(5,6) * t91 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t93 + (Ifges(7,4) + Ifges(6,5)) * t81 + (-Ifges(6,6) + Ifges(7,6)) * t80) * t93 + (t102 * (-mrSges(3,1) * t107 + mrSges(3,2) * t106) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t118 + mrSges(3,3)) * (t106 ^ 2 + t120) * qJ(2) + Ifges(3,2) * t120 / 0.2e1 + (Ifges(3,4) * t107 + Ifges(3,1) * t106 / 0.2e1) * t106) * qJD(1)) * qJD(1);
T  = t1;
