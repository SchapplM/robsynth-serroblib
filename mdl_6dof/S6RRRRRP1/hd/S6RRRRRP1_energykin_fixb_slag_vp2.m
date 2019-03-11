% Calculate kinetic energy for
% S6RRRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d3,d4,d5]';
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
% Datum: 2019-03-10 00:59
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRRRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRRRRP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-10 00:55:30
% EndTime: 2019-03-10 00:55:30
% DurationCPUTime: 0.62s
% Computational Cost: add. (798->106), mult. (1745->155), div. (0->0), fcn. (1288->8), ass. (0->38)
t118 = qJD(1) * (-pkin(8) - pkin(7));
t116 = pkin(7) * mrSges(3,3);
t106 = sin(qJ(5));
t110 = cos(qJ(5));
t105 = qJD(2) + qJD(3);
t104 = qJD(4) + t105;
t107 = sin(qJ(4));
t111 = cos(qJ(4));
t109 = sin(qJ(2));
t101 = qJD(2) * pkin(2) + t109 * t118;
t113 = cos(qJ(2));
t102 = t113 * t118;
t108 = sin(qJ(3));
t112 = cos(qJ(3));
t93 = t112 * t101 + t102 * t108;
t99 = (t108 * t113 + t109 * t112) * qJD(1);
t84 = pkin(3) * t105 - pkin(9) * t99 + t93;
t94 = t108 * t101 - t112 * t102;
t98 = (-t108 * t109 + t112 * t113) * qJD(1);
t87 = pkin(9) * t98 + t94;
t79 = t107 * t84 + t111 * t87;
t77 = pkin(10) * t104 + t79;
t91 = -t107 * t99 + t111 * t98;
t92 = t107 * t98 + t111 * t99;
t103 = (-pkin(2) * t113 - pkin(1)) * qJD(1);
t95 = -pkin(3) * t98 + t103;
t82 = -pkin(4) * t91 - pkin(10) * t92 + t95;
t73 = t106 * t82 + t110 * t77;
t72 = -t106 * t77 + t110 * t82;
t78 = -t107 * t87 + t111 * t84;
t76 = -pkin(4) * t104 - t78;
t90 = qJD(5) - t91;
t89 = t104 * t106 + t110 * t92;
t88 = t104 * t110 - t106 * t92;
t74 = -pkin(5) * t88 + qJD(6) + t76;
t71 = qJ(6) * t88 + t73;
t70 = pkin(5) * t90 - qJ(6) * t89 + t72;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t103 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t76 ^ 2) / 0.2e1 + (t103 * mrSges(4,2) - t93 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t95 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 + (-t103 * mrSges(4,1) + t94 * mrSges(4,3) + Ifges(4,4) * t99 + Ifges(4,2) * t98 / 0.2e1) * t98 + (-t95 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t92 + Ifges(5,2) * t91 / 0.2e1) * t91 + (t93 * mrSges(4,1) - t94 * mrSges(4,2) + Ifges(4,5) * t99 + Ifges(4,6) * t98 + Ifges(4,3) * t105 / 0.2e1) * t105 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,5) * t92 + Ifges(5,6) * t91 + Ifges(5,3) * t104 / 0.2e1) * t104 + (t72 * mrSges(6,1) + t70 * mrSges(7,1) - t73 * mrSges(6,2) - t71 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t90) * t90 + (t76 * mrSges(6,2) + t74 * mrSges(7,2) - t72 * mrSges(6,3) - t70 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t89 + (Ifges(6,5) + Ifges(7,5)) * t90) * t89 + (-t76 * mrSges(6,1) - t74 * mrSges(7,1) + t73 * mrSges(6,3) + t71 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t88 + (Ifges(6,6) + Ifges(7,6)) * t90 + (Ifges(6,4) + Ifges(7,4)) * t89) * t88 + ((Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t109 ^ 2 + t113 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t116 + Ifges(3,2) / 0.2e1) * t113) * t113 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t113 + (t116 + Ifges(3,1) / 0.2e1) * t109) * t109) * qJD(1) + ((-pkin(7) * mrSges(3,2) + Ifges(3,6)) * t113 + (-pkin(7) * mrSges(3,1) + Ifges(3,5)) * t109) * qJD(2)) * qJD(1);
T  = t1;
