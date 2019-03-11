% Calculate kinetic energy for
% S6RPRPRR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:50
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:47:57
% EndTime: 2019-03-09 03:47:58
% DurationCPUTime: 0.58s
% Computational Cost: add. (527->100), mult. (1259->146), div. (0->0), fcn. (886->8), ass. (0->42)
t113 = cos(pkin(10));
t129 = t113 ^ 2;
t128 = m(3) / 0.2e1;
t127 = cos(qJ(3));
t126 = pkin(7) + qJ(2);
t116 = sin(qJ(5));
t119 = cos(qJ(5));
t112 = sin(pkin(10));
t117 = sin(qJ(3));
t101 = (t127 * t112 + t113 * t117) * qJD(1);
t123 = t112 * qJD(1);
t102 = t126 * t123;
t124 = qJD(1) * t113;
t103 = t126 * t124;
t93 = -t127 * t102 - t117 * t103;
t122 = qJD(4) - t93;
t82 = -t101 * pkin(8) + (-pkin(3) - pkin(4)) * qJD(3) + t122;
t100 = t117 * t123 - t127 * t124;
t94 = -t117 * t102 + t127 * t103;
t92 = qJD(3) * qJ(4) + t94;
t84 = pkin(8) * t100 + t92;
t79 = t116 * t82 + t119 * t84;
t107 = -qJD(1) * pkin(1) + qJD(2);
t104 = -pkin(2) * t124 + t107;
t87 = t100 * pkin(3) - t101 * qJ(4) + t104;
t78 = -t116 * t84 + t119 * t82;
t89 = t100 * t119 - t101 * t116;
t80 = -pkin(4) * t100 - t87;
t118 = cos(qJ(6));
t115 = sin(qJ(6));
t110 = -qJD(3) + qJD(5);
t91 = -qJD(3) * pkin(3) + t122;
t90 = t100 * t116 + t101 * t119;
t88 = qJD(6) - t89;
t86 = t110 * t115 + t118 * t90;
t85 = t110 * t118 - t115 * t90;
t77 = pkin(9) * t110 + t79;
t76 = -pkin(5) * t110 - t78;
t75 = -pkin(5) * t89 - pkin(9) * t90 + t80;
t74 = t115 * t75 + t118 * t77;
t73 = -t115 * t77 + t118 * t75;
t1 = m(5) * (t87 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + t107 ^ 2 * t128 + m(4) * (t104 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + (t80 * mrSges(6,2) - t78 * mrSges(6,3) + Ifges(6,1) * t90 / 0.2e1) * t90 + (t73 * mrSges(7,1) - t74 * mrSges(7,2) + Ifges(7,3) * t88 / 0.2e1) * t88 + (-t80 * mrSges(6,1) + t79 * mrSges(6,3) + Ifges(6,4) * t90 + Ifges(6,2) * t89 / 0.2e1) * t89 + (t76 * mrSges(7,2) - t73 * mrSges(7,3) + Ifges(7,5) * t88 + Ifges(7,1) * t86 / 0.2e1) * t86 + (-t76 * mrSges(7,1) + t74 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,6) * t88 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t78 * mrSges(6,1) - t79 * mrSges(6,2) + Ifges(6,5) * t90 + Ifges(6,6) * t89 + Ifges(6,3) * t110 / 0.2e1) * t110 + (t104 * mrSges(4,2) + t91 * mrSges(5,2) - t93 * mrSges(4,3) - t87 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * t101) * t101 + (t104 * mrSges(4,1) + t87 * mrSges(5,1) - t92 * mrSges(5,2) - t94 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t100 + (-Ifges(4,4) + Ifges(5,5)) * t101) * t100 + (t93 * mrSges(4,1) - t91 * mrSges(5,1) - t94 * mrSges(4,2) + t92 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (Ifges(5,4) + Ifges(4,5)) * t101 + (-Ifges(4,6) + Ifges(5,6)) * t100) * qJD(3) + (t107 * (-mrSges(3,1) * t113 + mrSges(3,2) * t112) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t128 + mrSges(3,3)) * (t112 ^ 2 + t129) * qJ(2) + Ifges(3,2) * t129 / 0.2e1 + (Ifges(3,4) * t113 + Ifges(3,1) * t112 / 0.2e1) * t112) * qJD(1)) * qJD(1);
T  = t1;
