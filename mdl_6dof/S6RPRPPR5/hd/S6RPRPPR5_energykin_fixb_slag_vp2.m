% Calculate kinetic energy for
% S6RPRPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2,theta5]';
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
% Datum: 2019-03-09 02:52
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:50:03
% EndTime: 2019-03-09 02:50:03
% DurationCPUTime: 0.54s
% Computational Cost: add. (531->101), mult. (1269->145), div. (0->0), fcn. (888->8), ass. (0->43)
t109 = cos(pkin(9));
t124 = t109 ^ 2;
t123 = m(3) / 0.2e1;
t122 = cos(qJ(3));
t121 = pkin(3) + qJ(5);
t120 = pkin(7) + qJ(2);
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t100 = qJD(2) + (-pkin(2) * t109 - pkin(1)) * qJD(1);
t107 = sin(pkin(9));
t111 = sin(qJ(3));
t97 = (t107 * t122 + t109 * t111) * qJD(1);
t115 = -qJ(4) * t97 + t100;
t117 = qJD(1) * t109;
t118 = qJD(1) * t107;
t96 = t111 * t118 - t117 * t122;
t78 = t121 * t96 + t115;
t98 = t120 * t118;
t99 = t120 * t117;
t88 = -t111 * t99 - t122 * t98;
t116 = qJD(4) - t88;
t81 = t97 * pkin(4) - qJD(3) * t121 + t116;
t75 = t106 * t81 + t108 * t78;
t89 = -t111 * t98 + t122 * t99;
t87 = -qJD(3) * qJ(4) - t89;
t74 = -t106 * t78 + t108 * t81;
t84 = -pkin(4) * t96 + qJD(5) - t87;
t112 = cos(qJ(6));
t110 = sin(qJ(6));
t102 = -qJD(1) * pkin(1) + qJD(2);
t92 = qJD(6) + t97;
t91 = qJD(3) * t108 + t106 * t96;
t90 = -qJD(3) * t106 + t108 * t96;
t86 = -qJD(3) * pkin(3) + t116;
t85 = pkin(3) * t96 + t115;
t83 = t110 * t90 + t112 * t91;
t82 = -t110 * t91 + t112 * t90;
t76 = -pkin(5) * t90 + t84;
t73 = pkin(8) * t90 + t75;
t72 = pkin(5) * t97 - pkin(8) * t91 + t74;
t71 = t110 * t72 + t112 * t73;
t70 = -t110 * t73 + t112 * t72;
t1 = m(4) * (t100 ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + t102 ^ 2 * t123 + m(6) * (t74 ^ 2 + t75 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t76 ^ 2) / 0.2e1 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t92 / 0.2e1) * t92 + (t84 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,1) * t91 / 0.2e1) * t91 + (-t84 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,2) * t90 / 0.2e1) * t90 + (t76 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t92 + Ifges(7,1) * t83 / 0.2e1) * t83 + (-t76 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t83 + Ifges(7,6) * t92 + Ifges(7,2) * t82 / 0.2e1) * t82 + (t100 * mrSges(4,1) + t87 * mrSges(5,1) - t85 * mrSges(5,2) - t89 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t96) * t96 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + t86 * mrSges(5,2) - t87 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (Ifges(5,5) - Ifges(4,6)) * t96) * qJD(3) + (t86 * mrSges(5,1) + t74 * mrSges(6,1) + t100 * mrSges(4,2) - t75 * mrSges(6,2) - t88 * mrSges(4,3) - t85 * mrSges(5,3) + Ifges(6,5) * t91 + Ifges(6,6) * t90 + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t97 + (-Ifges(4,4) - Ifges(5,6)) * t96 + (-Ifges(5,4) + Ifges(4,5)) * qJD(3)) * t97 + (t102 * (-mrSges(3,1) * t109 + mrSges(3,2) * t107) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t123 + mrSges(3,3)) * (t107 ^ 2 + t124) * qJ(2) + Ifges(3,2) * t124 / 0.2e1 + (Ifges(3,4) * t109 + Ifges(3,1) * t107 / 0.2e1) * t107) * qJD(1)) * qJD(1);
T  = t1;
