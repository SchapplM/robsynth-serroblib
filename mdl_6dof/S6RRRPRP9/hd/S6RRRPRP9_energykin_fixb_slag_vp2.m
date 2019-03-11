% Calculate kinetic energy for
% S6RRRPRP9
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d5]';
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
% Datum: 2019-03-09 17:28
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRPRP9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRPRP9_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 17:21:47
% EndTime: 2019-03-09 17:21:47
% DurationCPUTime: 0.48s
% Computational Cost: add. (494->104), mult. (971->139), div. (0->0), fcn. (608->6), ass. (0->33)
t112 = pkin(7) * mrSges(3,3);
t111 = cos(qJ(3));
t110 = cos(qJ(5));
t101 = sin(qJ(5));
t102 = sin(qJ(3));
t103 = sin(qJ(2));
t104 = cos(qJ(2));
t86 = (-pkin(2) * t104 - pkin(8) * t103 - pkin(1)) * qJD(1);
t109 = qJD(1) * t104;
t94 = pkin(7) * t109 + qJD(2) * pkin(8);
t82 = -t102 * t94 + t111 * t86;
t107 = qJD(4) - t82;
t108 = t103 * qJD(1);
t90 = t102 * qJD(2) + t108 * t111;
t97 = qJD(3) - t109;
t73 = -t90 * pkin(9) + (-pkin(3) - pkin(4)) * t97 + t107;
t83 = t102 * t86 + t111 * t94;
t78 = t97 * qJ(4) + t83;
t89 = -qJD(2) * t111 + t102 * t108;
t75 = pkin(9) * t89 + t78;
t70 = t101 * t73 + t110 * t75;
t93 = -qJD(2) * pkin(2) + pkin(7) * t108;
t79 = t89 * pkin(3) - t90 * qJ(4) + t93;
t76 = -pkin(4) * t89 - t79;
t69 = -t101 * t75 + t110 * t73;
t96 = qJD(5) - t97;
t81 = t101 * t89 + t110 * t90;
t80 = t101 * t90 - t110 * t89;
t77 = -t97 * pkin(3) + t107;
t71 = pkin(5) * t80 - qJ(6) * t81 + t76;
t68 = qJ(6) * t96 + t70;
t67 = -t96 * pkin(5) + qJD(6) - t69;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t82 ^ 2 + t83 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t76 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + (t82 * mrSges(4,1) - t77 * mrSges(5,1) - t83 * mrSges(4,2) + t78 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t97) * t97 + (t69 * mrSges(6,1) - t67 * mrSges(7,1) - t70 * mrSges(6,2) + t68 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t96) * t96 + (t93 * mrSges(4,2) + t77 * mrSges(5,2) - t82 * mrSges(4,3) - t79 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t90 + (Ifges(5,4) + Ifges(4,5)) * t97) * t90 + (t76 * mrSges(6,2) + t67 * mrSges(7,2) - t69 * mrSges(6,3) - t71 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t81 + (Ifges(7,4) + Ifges(6,5)) * t96) * t81 + (t93 * mrSges(4,1) + t79 * mrSges(5,1) - t78 * mrSges(5,2) - t83 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t89 + (-Ifges(4,6) + Ifges(5,6)) * t97 + (-Ifges(4,4) + Ifges(5,5)) * t90) * t89 + (t76 * mrSges(6,1) + t71 * mrSges(7,1) - t68 * mrSges(7,2) - t70 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t96 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t103 ^ 2 + t104 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t112 + Ifges(3,2) / 0.2e1) * t104) * t104 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t104 + (t112 + Ifges(3,1) / 0.2e1) * t103) * t103) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t104 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t103) * qJD(2)) * qJD(1);
T  = t1;
