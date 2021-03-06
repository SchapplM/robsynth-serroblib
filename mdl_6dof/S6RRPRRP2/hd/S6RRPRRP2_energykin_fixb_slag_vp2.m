% Calculate kinetic energy for
% S6RRPRRP2
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
% Datum: 2019-03-09 11:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPRRP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 11:43:00
% EndTime: 2019-03-09 11:43:01
% DurationCPUTime: 0.60s
% Computational Cost: add. (738->106), mult. (1741->154), div. (0->0), fcn. (1284->8), ass. (0->37)
t117 = qJD(1) * (pkin(7) + qJ(3));
t115 = pkin(7) * mrSges(3,3);
t114 = cos(qJ(5));
t106 = sin(qJ(5));
t103 = qJD(2) + qJD(4);
t107 = sin(qJ(4));
t109 = cos(qJ(4));
t110 = cos(qJ(2));
t100 = t110 * t117;
t104 = sin(pkin(10));
t105 = cos(pkin(10));
t108 = sin(qJ(2));
t99 = qJD(2) * pkin(2) - t108 * t117;
t91 = -t100 * t104 + t105 * t99;
t97 = (t104 * t110 + t105 * t108) * qJD(1);
t84 = qJD(2) * pkin(3) - pkin(8) * t97 + t91;
t92 = t105 * t100 + t104 * t99;
t96 = (-t104 * t108 + t105 * t110) * qJD(1);
t85 = pkin(8) * t96 + t92;
t78 = t107 * t84 + t109 * t85;
t76 = pkin(9) * t103 + t78;
t89 = -t107 * t97 + t109 * t96;
t90 = t107 * t96 + t109 * t97;
t101 = qJD(3) + (-pkin(2) * t110 - pkin(1)) * qJD(1);
t93 = -pkin(3) * t96 + t101;
t80 = -pkin(4) * t89 - pkin(9) * t90 + t93;
t72 = t106 * t80 + t114 * t76;
t77 = -t107 * t85 + t109 * t84;
t75 = -pkin(4) * t103 - t77;
t71 = -t106 * t76 + t114 * t80;
t88 = qJD(5) - t89;
t87 = t106 * t103 + t114 * t90;
t86 = -t114 * t103 + t106 * t90;
t73 = pkin(5) * t86 - qJ(6) * t87 + t75;
t70 = qJ(6) * t88 + t72;
t69 = -t88 * pkin(5) + qJD(6) - t71;
t1 = m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t91 ^ 2 + t92 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t93 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t91 * mrSges(4,3) + Ifges(4,1) * t97 / 0.2e1) * t97 + (t93 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t90 / 0.2e1) * t90 + (-t101 * mrSges(4,1) + t92 * mrSges(4,3) + Ifges(4,4) * t97 + Ifges(4,2) * t96 / 0.2e1) * t96 + (-t93 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t90 + Ifges(5,2) * t89 / 0.2e1) * t89 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + Ifges(5,5) * t90 + Ifges(5,6) * t89 + Ifges(5,3) * t103 / 0.2e1) * t103 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,5) * t97 + Ifges(4,6) * t96 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(2) + (Ifges(3,5) * t108 + Ifges(3,6) * t110 + (-mrSges(3,1) * t108 - mrSges(3,2) * t110) * pkin(7)) * qJD(1)) * qJD(2) + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t88) * t88 + (t75 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t87 + (Ifges(7,4) + Ifges(6,5)) * t88) * t87 + (t75 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t86 + (-Ifges(6,6) + Ifges(7,6)) * t88 + (-Ifges(6,4) + Ifges(7,5)) * t87) * t86 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t108 ^ 2 + t110 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (pkin(1) * mrSges(3,1) + (Ifges(3,2) / 0.2e1 + t115) * t110) * t110 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t110 + (Ifges(3,1) / 0.2e1 + t115) * t108) * t108) * qJD(1) ^ 2;
T  = t1;
