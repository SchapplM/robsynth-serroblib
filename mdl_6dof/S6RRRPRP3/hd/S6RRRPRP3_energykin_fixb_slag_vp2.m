% Calculate kinetic energy for
% S6RRRPRP3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-03-09 16:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 16:39:02
% EndTime: 2019-03-09 16:39:03
% DurationCPUTime: 0.55s
% Computational Cost: add. (792->106), mult. (1641->153), div. (0->0), fcn. (1186->8), ass. (0->39)
t119 = -pkin(8) - pkin(7);
t118 = pkin(7) * mrSges(3,3);
t117 = cos(qJ(5));
t108 = sin(qJ(5));
t106 = sin(pkin(10));
t107 = cos(pkin(10));
t112 = cos(qJ(2));
t103 = (-pkin(2) * t112 - pkin(1)) * qJD(1);
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t115 = qJD(1) * t112;
t110 = sin(qJ(2));
t116 = qJD(1) * t110;
t97 = t109 * t116 - t111 * t115;
t98 = (t109 * t112 + t110 * t111) * qJD(1);
t86 = pkin(3) * t97 - qJ(4) * t98 + t103;
t105 = qJD(2) + qJD(3);
t101 = qJD(2) * pkin(2) + t119 * t116;
t102 = t119 * t115;
t91 = t109 * t101 - t111 * t102;
t89 = qJ(4) * t105 + t91;
t79 = -t106 * t89 + t107 * t86;
t94 = t105 * t106 + t107 * t98;
t76 = pkin(4) * t97 - pkin(9) * t94 + t79;
t80 = t106 * t86 + t107 * t89;
t93 = t105 * t107 - t106 * t98;
t78 = pkin(9) * t93 + t80;
t73 = t108 * t76 + t117 * t78;
t90 = t101 * t111 + t109 * t102;
t72 = -t108 * t78 + t117 * t76;
t88 = -pkin(3) * t105 + qJD(4) - t90;
t81 = -pkin(4) * t93 + t88;
t96 = qJD(5) + t97;
t83 = t108 * t93 + t117 * t94;
t82 = t108 * t94 - t117 * t93;
t74 = pkin(5) * t82 - qJ(6) * t83 + t81;
t71 = qJ(6) * t96 + t73;
t70 = -t96 * pkin(5) + qJD(6) - t72;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t103 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t79 ^ 2 + t80 ^ 2 + t88 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t81 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) / 0.2e1 + (t103 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t88 * mrSges(5,2) - t79 * mrSges(5,3) + Ifges(5,1) * t94 / 0.2e1) * t94 + (-t88 * mrSges(5,1) + t80 * mrSges(5,3) + Ifges(5,4) * t94 + Ifges(5,2) * t93 / 0.2e1) * t93 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,5) * t98 + Ifges(4,3) * t105 / 0.2e1) * t105 + (t72 * mrSges(6,1) - t70 * mrSges(7,1) - t73 * mrSges(6,2) + t71 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t96) * t96 + (t81 * mrSges(6,2) + t70 * mrSges(7,2) - t72 * mrSges(6,3) - t74 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t83 + (Ifges(7,4) + Ifges(6,5)) * t96) * t83 + (t103 * mrSges(4,1) + t79 * mrSges(5,1) - t80 * mrSges(5,2) - t91 * mrSges(4,3) - Ifges(4,4) * t98 + Ifges(5,5) * t94 - Ifges(4,6) * t105 + Ifges(5,6) * t93 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t97) * t97 + (t81 * mrSges(6,1) + t74 * mrSges(7,1) - t71 * mrSges(7,2) - t73 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t82 + (-Ifges(6,6) + Ifges(7,6)) * t96 + (-Ifges(6,4) + Ifges(7,5)) * t83) * t82 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t110 ^ 2 + t112 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t118) * t112) * t112 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t112 + (Ifges(3,1) / 0.2e1 + t118) * t110) * t110) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t112 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t110) * qJD(2)) * qJD(1);
T  = t1;
