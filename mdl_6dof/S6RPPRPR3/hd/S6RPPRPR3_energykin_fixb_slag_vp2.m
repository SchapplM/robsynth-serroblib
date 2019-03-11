% Calculate kinetic energy for
% S6RPPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta5]';
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
% Datum: 2019-03-09 01:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:43:59
% EndTime: 2019-03-09 01:43:59
% DurationCPUTime: 0.36s
% Computational Cost: add. (302->82), mult. (598->126), div. (0->0), fcn. (330->8), ass. (0->34)
t100 = sin(qJ(4));
t102 = cos(qJ(4));
t95 = sin(pkin(10));
t97 = cos(pkin(10));
t87 = (t100 * t97 + t102 * t95) * qJD(1);
t96 = sin(pkin(9));
t91 = (pkin(1) * t96 + qJ(3)) * qJD(1);
t111 = t91 ^ 2;
t110 = m(3) / 0.2e1;
t98 = cos(pkin(9));
t108 = -pkin(1) * t98 - pkin(2);
t89 = qJD(3) + (-pkin(7) + t108) * qJD(1);
t81 = -qJD(2) * t100 + t102 * t89;
t77 = -qJ(5) * qJD(1) * t102 + qJD(4) * pkin(4) + t81;
t109 = qJD(1) * t100;
t82 = qJD(2) * t102 + t100 * t89;
t78 = -qJ(5) * t109 + t82;
t73 = t77 * t95 + t78 * t97;
t72 = t77 * t97 - t78 * t95;
t86 = pkin(4) * t109 + qJD(5) + t91;
t103 = qJD(2) ^ 2;
t101 = cos(qJ(6));
t99 = sin(qJ(6));
t90 = qJD(1) * t108 + qJD(3);
t88 = (-t100 * t95 + t102 * t97) * qJD(1);
t83 = qJD(6) + t87;
t80 = qJD(4) * t99 + t101 * t88;
t79 = qJD(4) * t101 - t88 * t99;
t74 = t87 * pkin(5) - t88 * pkin(8) + t86;
t71 = qJD(4) * pkin(8) + t73;
t70 = -qJD(4) * pkin(5) - t72;
t69 = t101 * t71 + t74 * t99;
t68 = t101 * t74 - t71 * t99;
t1 = m(4) * (t90 ^ 2 + t103 + t111) / 0.2e1 + t103 * t110 + m(5) * (t81 ^ 2 + t82 ^ 2 + t111) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (t86 * mrSges(6,2) - t72 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t83 / 0.2e1) * t83 - (-t86 * mrSges(6,1) + t73 * mrSges(6,3) + Ifges(6,4) * t88 - Ifges(6,2) * t87 / 0.2e1) * t87 + (t70 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t83 + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t70 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,6) * t83 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t81 * mrSges(5,1) + t72 * mrSges(6,1) - t82 * mrSges(5,2) - t73 * mrSges(6,2) + Ifges(6,5) * t88 - Ifges(6,6) * t87 + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(4)) * qJD(4) + (t90 * mrSges(4,2) + (mrSges(5,1) * t100 + mrSges(5,2) * t102 + mrSges(4,3)) * t91 + (-t100 * t82 - t102 * t81) * mrSges(5,3) + qJD(4) * (Ifges(5,5) * t102 - Ifges(5,6) * t100) + (Ifges(2,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t98 * mrSges(3,1) - t96 * mrSges(3,2) + (t96 ^ 2 + t98 ^ 2) * t110 * pkin(1)) * pkin(1) + Ifges(5,1) * t102 ^ 2 / 0.2e1 + (-Ifges(5,4) * t102 + Ifges(5,2) * t100 / 0.2e1) * t100) * qJD(1)) * qJD(1);
T  = t1;
