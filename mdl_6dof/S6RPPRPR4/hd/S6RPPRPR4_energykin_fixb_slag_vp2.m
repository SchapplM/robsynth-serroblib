% Calculate kinetic energy for
% S6RPPRPR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:46:11
% EndTime: 2019-03-09 01:46:12
% DurationCPUTime: 0.36s
% Computational Cost: add. (357->84), mult. (670->129), div. (0->0), fcn. (362->8), ass. (0->36)
t100 = sin(pkin(10));
t102 = cos(pkin(10));
t106 = sin(qJ(4));
t108 = cos(qJ(4));
t90 = (t100 * t106 - t102 * t108) * qJD(1);
t113 = m(3) / 0.2e1;
t101 = sin(pkin(9));
t103 = cos(pkin(9));
t112 = qJ(2) * qJD(1);
t93 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t88 = t101 * t93 + t103 * t112;
t86 = -qJD(1) * pkin(7) + t88;
t99 = t108 * qJD(3);
t77 = qJD(4) * pkin(4) + t99 + (qJ(5) * qJD(1) - t86) * t106;
t111 = qJD(1) * t108;
t80 = t106 * qJD(3) + t108 * t86;
t78 = -qJ(5) * t111 + t80;
t73 = t100 * t77 + t102 * t78;
t87 = -t101 * t112 + t103 * t93;
t85 = qJD(1) * pkin(3) - t87;
t72 = -t100 * t78 + t102 * t77;
t81 = pkin(4) * t111 + qJD(5) + t85;
t107 = cos(qJ(6));
t105 = sin(qJ(6));
t97 = -qJD(1) * pkin(1) + qJD(2);
t91 = (-t100 * t108 - t102 * t106) * qJD(1);
t89 = qJD(6) - t90;
t84 = t105 * qJD(4) + t107 * t91;
t83 = t107 * qJD(4) - t105 * t91;
t79 = -t106 * t86 + t99;
t74 = -t90 * pkin(5) - t91 * pkin(8) + t81;
t71 = qJD(4) * pkin(8) + t73;
t70 = -qJD(4) * pkin(5) - t72;
t69 = t105 * t74 + t107 * t71;
t68 = -t105 * t71 + t107 * t74;
t1 = t97 ^ 2 * t113 + m(5) * (t79 ^ 2 + t80 ^ 2 + t85 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t87 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (t81 * mrSges(6,2) - t72 * mrSges(6,3) + Ifges(6,1) * t91 / 0.2e1) * t91 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t89 / 0.2e1) * t89 + (-t81 * mrSges(6,1) + t73 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t90 / 0.2e1) * t90 + (t70 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t89 + Ifges(7,1) * t84 / 0.2e1) * t84 + (-t70 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t84 + Ifges(7,6) * t89 + Ifges(7,2) * t83 / 0.2e1) * t83 + (t79 * mrSges(5,1) + t72 * mrSges(6,1) - t80 * mrSges(5,2) - t73 * mrSges(6,2) + Ifges(6,5) * t91 + Ifges(6,6) * t90 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(4)) * qJD(4) + (t85 * (mrSges(5,1) * t108 - mrSges(5,2) * t106) - t87 * mrSges(4,1) + t88 * mrSges(4,2) - t97 * mrSges(3,1) + (t79 * t106 - t80 * t108) * mrSges(5,3) + qJD(4) * (-Ifges(5,5) * t106 - Ifges(5,6) * t108) + (Ifges(3,2) / 0.2e1 + Ifges(2,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t113 + mrSges(3,3)) * qJ(2) + t108 ^ 2 * Ifges(5,2) / 0.2e1 + (Ifges(5,4) * t108 + Ifges(5,1) * t106 / 0.2e1) * t106) * qJD(1)) * qJD(1);
T  = t1;
