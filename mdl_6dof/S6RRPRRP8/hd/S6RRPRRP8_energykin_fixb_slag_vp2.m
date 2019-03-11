% Calculate kinetic energy for
% S6RRPRRP8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-03-09 12:26
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 12:21:52
% EndTime: 2019-03-09 12:21:53
% DurationCPUTime: 0.59s
% Computational Cost: add. (810->107), mult. (1787->155), div. (0->0), fcn. (1290->8), ass. (0->38)
t117 = pkin(7) * mrSges(3,3);
t116 = cos(qJ(5));
t107 = sin(qJ(5));
t111 = cos(qJ(2));
t114 = t111 * qJD(1);
t103 = qJD(4) - t114;
t108 = sin(qJ(4));
t110 = cos(qJ(4));
t101 = pkin(7) * t114 + qJD(2) * qJ(3);
t105 = sin(pkin(10));
t106 = cos(pkin(10));
t109 = sin(qJ(2));
t96 = (-pkin(2) * t111 - qJ(3) * t109 - pkin(1)) * qJD(1);
t91 = -t101 * t105 + t106 * t96;
t115 = t109 * qJD(1);
t98 = qJD(2) * t105 + t106 * t115;
t85 = -pkin(3) * t114 - pkin(8) * t98 + t91;
t92 = t106 * t101 + t105 * t96;
t97 = qJD(2) * t106 - t105 * t115;
t87 = pkin(8) * t97 + t92;
t78 = -t108 * t87 + t110 * t85;
t90 = t108 * t97 + t110 * t98;
t75 = pkin(4) * t103 - pkin(9) * t90 + t78;
t79 = t108 * t85 + t110 * t87;
t89 = -t108 * t98 + t110 * t97;
t77 = pkin(9) * t89 + t79;
t72 = t107 * t75 + t116 * t77;
t100 = -qJD(2) * pkin(2) + pkin(7) * t115 + qJD(3);
t71 = -t107 * t77 + t116 * t75;
t93 = -pkin(3) * t97 + t100;
t82 = -pkin(4) * t89 + t93;
t102 = qJD(5) + t103;
t81 = t107 * t89 + t116 * t90;
t80 = t107 * t90 - t116 * t89;
t73 = pkin(5) * t80 - qJ(6) * t81 + t82;
t70 = qJ(6) * t102 + t72;
t69 = -t102 * pkin(5) + qJD(6) - t71;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t100 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + (t100 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t93 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t90 / 0.2e1) * t90 + (-t93 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t90 + Ifges(5,2) * t89 / 0.2e1) * t89 + (-Ifges(4,6) * t114 - t100 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,4) * t98 + Ifges(4,2) * t97 / 0.2e1) * t97 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,5) * t90 + Ifges(5,6) * t89 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t82 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81) * t81 + (t82 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,3) / 0.2e1) * t80 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t102 + (Ifges(7,4) + Ifges(6,5)) * t81 + (-Ifges(6,6) + Ifges(7,6)) * t80) * t102 + ((-pkin(7) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t109 + (-t91 * mrSges(4,1) + t92 * mrSges(4,2) - Ifges(4,5) * t98 + (-pkin(7) * mrSges(3,2) + Ifges(3,6)) * qJD(2)) * t111 + ((-pkin(1) * mrSges(3,2) + (t117 + Ifges(3,1) / 0.2e1) * t109) * t109 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t109 + (t117 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t111) * t111 + Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t109 ^ 2 + t111 ^ 2) * pkin(7) ^ 2) / 0.2e1) * qJD(1)) * qJD(1);
T  = t1;
