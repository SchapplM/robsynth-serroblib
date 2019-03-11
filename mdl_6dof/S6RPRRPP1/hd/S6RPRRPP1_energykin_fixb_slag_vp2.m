% Calculate kinetic energy for
% S6RPRRPP1
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
% Datum: 2019-03-09 04:30
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:27:49
% EndTime: 2019-03-09 04:27:50
% DurationCPUTime: 0.44s
% Computational Cost: add. (442->93), mult. (919->137), div. (0->0), fcn. (574->8), ass. (0->34)
t110 = m(3) / 0.2e1;
t109 = cos(pkin(10));
t101 = sin(qJ(4));
t103 = cos(qJ(4));
t102 = sin(qJ(3));
t104 = cos(qJ(3));
t99 = sin(pkin(9));
t94 = (pkin(1) * t99 + pkin(7)) * qJD(1);
t88 = t102 * qJD(2) + t104 * t94;
t85 = qJD(3) * pkin(8) + t88;
t100 = cos(pkin(9));
t107 = -pkin(1) * t100 - pkin(2);
t86 = (-pkin(3) * t104 - pkin(8) * t102 + t107) * qJD(1);
t76 = -t101 * t85 + t103 * t86;
t108 = qJD(1) * t102;
t93 = qJD(3) * t101 + t103 * t108;
t96 = -qJD(1) * t104 + qJD(4);
t73 = pkin(4) * t96 - qJ(5) * t93 + t76;
t77 = t101 * t86 + t103 * t85;
t92 = qJD(3) * t103 - t101 * t108;
t75 = qJ(5) * t92 + t77;
t98 = sin(pkin(10));
t70 = t109 * t75 + t98 * t73;
t87 = qJD(2) * t104 - t102 * t94;
t84 = -qJD(3) * pkin(3) - t87;
t69 = t109 * t73 - t98 * t75;
t78 = -pkin(4) * t92 + qJD(5) + t84;
t95 = t107 * qJD(1);
t80 = t109 * t93 + t98 * t92;
t79 = -t109 * t92 + t93 * t98;
t71 = pkin(5) * t79 - qJ(6) * t80 + t78;
t68 = qJ(6) * t96 + t70;
t67 = -t96 * pkin(5) + qJD(6) - t69;
t1 = qJD(2) ^ 2 * t110 + m(4) * (t87 ^ 2 + t88 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t78 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + (t84 * mrSges(5,2) - t76 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (-t84 * mrSges(5,1) + t77 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (t87 * mrSges(4,1) - t88 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t78 * mrSges(6,2) + t67 * mrSges(7,2) - t69 * mrSges(6,3) - t71 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t80) * t80 + (t78 * mrSges(6,1) + t71 * mrSges(7,1) - t68 * mrSges(7,2) - t70 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t79 + (-Ifges(6,4) + Ifges(7,5)) * t80) * t79 + (t76 * mrSges(5,1) + t69 * mrSges(6,1) - t67 * mrSges(7,1) - t77 * mrSges(5,2) - t70 * mrSges(6,2) + t68 * mrSges(7,3) + Ifges(5,5) * t93 + Ifges(5,6) * t92 + (Ifges(5,3) / 0.2e1 + Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t96 + (Ifges(7,4) + Ifges(6,5)) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t79) * t96 + (t95 * (-mrSges(4,1) * t104 + mrSges(4,2) * t102) + (-t87 * t102 + t88 * t104) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t102 + Ifges(4,6) * t104) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t100 * mrSges(3,1) - t99 * mrSges(3,2) + (t100 ^ 2 + t99 ^ 2) * t110 * pkin(1)) * pkin(1) + Ifges(4,2) * t104 ^ 2 / 0.2e1 + (Ifges(4,4) * t104 + Ifges(4,1) * t102 / 0.2e1) * t102) * qJD(1)) * qJD(1);
T  = t1;
