% Calculate kinetic energy for
% S6RPRPPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta4]';
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
% Datum: 2019-03-09 02:43
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:41:10
% EndTime: 2019-03-09 02:41:11
% DurationCPUTime: 0.44s
% Computational Cost: add. (352->95), mult. (783->139), div. (0->0), fcn. (474->8), ass. (0->38)
t112 = m(3) / 0.2e1;
t111 = pkin(4) + pkin(8);
t110 = cos(pkin(10));
t101 = sin(qJ(3));
t98 = sin(pkin(9));
t92 = (pkin(1) * t98 + pkin(7)) * qJD(1);
t103 = cos(qJ(3));
t96 = t103 * qJD(2);
t81 = qJD(3) * pkin(3) + t96 + (-qJ(4) * qJD(1) - t92) * t101;
t109 = qJD(1) * t103;
t86 = t101 * qJD(2) + t103 * t92;
t84 = qJ(4) * t109 + t86;
t97 = sin(pkin(10));
t76 = t110 * t84 + t97 * t81;
t99 = cos(pkin(9));
t108 = -pkin(1) * t99 - pkin(2);
t74 = -qJD(3) * qJ(5) - t76;
t75 = t110 * t81 - t97 * t84;
t107 = qJD(5) - t75;
t90 = (t110 * t101 + t103 * t97) * qJD(1);
t88 = qJD(4) + (-pkin(3) * t103 + t108) * qJD(1);
t106 = -qJ(5) * t90 + t88;
t102 = cos(qJ(6));
t100 = sin(qJ(6));
t93 = t108 * qJD(1);
t89 = qJD(1) * t101 * t97 - t110 * t109;
t87 = qJD(6) + t90;
t85 = -t101 * t92 + t96;
t83 = qJD(3) * t102 + t100 * t89;
t82 = -qJD(3) * t100 + t102 * t89;
t77 = pkin(4) * t89 + t106;
t73 = -qJD(3) * pkin(4) + t107;
t72 = t111 * t89 + t106;
t71 = -pkin(5) * t89 - t74;
t70 = t90 * pkin(5) - t111 * qJD(3) + t107;
t69 = t100 * t70 + t102 * t72;
t68 = -t100 * t72 + t102 * t70;
t1 = qJD(2) ^ 2 * t112 + m(4) * (t85 ^ 2 + t86 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t88 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + (t68 * mrSges(7,1) - t69 * mrSges(7,2) + Ifges(7,3) * t87 / 0.2e1) * t87 + (t71 * mrSges(7,2) - t68 * mrSges(7,3) + Ifges(7,5) * t87 + Ifges(7,1) * t83 / 0.2e1) * t83 + (-t71 * mrSges(7,1) + t69 * mrSges(7,3) + Ifges(7,4) * t83 + Ifges(7,6) * t87 + Ifges(7,2) * t82 / 0.2e1) * t82 + (t73 * mrSges(6,1) + t88 * mrSges(5,2) - t75 * mrSges(5,3) - t77 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t90) * t90 + (t88 * mrSges(5,1) + t74 * mrSges(6,1) - t77 * mrSges(6,2) - t76 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t89 + (-Ifges(5,4) - Ifges(6,6)) * t90) * t89 + (t85 * mrSges(4,1) + t75 * mrSges(5,1) - t86 * mrSges(4,2) - t76 * mrSges(5,2) + t73 * mrSges(6,2) - t74 * mrSges(6,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(3) + (-Ifges(6,4) + Ifges(5,5)) * t90 + (Ifges(6,5) - Ifges(5,6)) * t89) * qJD(3) + (t93 * (-mrSges(4,1) * t103 + mrSges(4,2) * t101) + (-t85 * t101 + t86 * t103) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t101 + Ifges(4,6) * t103) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t99 * mrSges(3,1) - t98 * mrSges(3,2) + (t98 ^ 2 + t99 ^ 2) * t112 * pkin(1)) * pkin(1) + Ifges(4,2) * t103 ^ 2 / 0.2e1 + (Ifges(4,4) * t103 + Ifges(4,1) * t101 / 0.2e1) * t101) * qJD(1)) * qJD(1);
T  = t1;
