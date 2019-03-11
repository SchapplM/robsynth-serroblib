% Calculate kinetic energy for
% S6RPRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:08
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:04:21
% EndTime: 2019-03-09 05:04:22
% DurationCPUTime: 0.50s
% Computational Cost: add. (380->94), mult. (755->139), div. (0->0), fcn. (440->8), ass. (0->38)
t121 = m(3) / 0.2e1;
t120 = -pkin(4) - pkin(5);
t119 = cos(qJ(4));
t108 = sin(qJ(4));
t109 = sin(qJ(3));
t111 = cos(qJ(3));
t105 = sin(pkin(10));
t97 = (pkin(1) * t105 + pkin(7)) * qJD(1);
t92 = t109 * qJD(2) + t111 * t97;
t89 = qJD(3) * pkin(8) + t92;
t106 = cos(pkin(10));
t117 = -pkin(1) * t106 - pkin(2);
t90 = (-pkin(3) * t111 - pkin(8) * t109 + t117) * qJD(1);
t82 = t108 * t90 + t119 * t89;
t91 = t111 * qJD(2) - t109 * t97;
t118 = qJD(1) * t109;
t101 = -qJD(1) * t111 + qJD(4);
t79 = t101 * qJ(5) + t82;
t116 = qJD(3) * pkin(3) + t91;
t81 = -t108 * t89 + t119 * t90;
t115 = qJD(5) - t81;
t96 = t108 * qJD(3) + t119 * t118;
t114 = qJ(5) * t96 + t116;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t100 = qJD(6) - t101;
t98 = t117 * qJD(1);
t95 = -t119 * qJD(3) + t108 * t118;
t84 = t107 * t95 + t110 * t96;
t83 = -t107 * t96 + t110 * t95;
t80 = pkin(4) * t95 - t114;
t78 = -t101 * pkin(4) + t115;
t77 = t120 * t95 + t114;
t76 = pkin(9) * t95 + t79;
t75 = -t96 * pkin(9) + t120 * t101 + t115;
t74 = t107 * t75 + t110 * t76;
t73 = -t107 * t76 + t110 * t75;
t1 = qJD(2) ^ 2 * t121 + m(4) * (t91 ^ 2 + t92 ^ 2 + t98 ^ 2) / 0.2e1 + m(5) * (t116 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t73 ^ 2 + t74 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t78 ^ 2 + t79 ^ 2 + t80 ^ 2) / 0.2e1 + (t77 * mrSges(7,2) - t73 * mrSges(7,3) + Ifges(7,1) * t84 / 0.2e1) * t84 + (-t77 * mrSges(7,1) + t74 * mrSges(7,3) + Ifges(7,4) * t84 + Ifges(7,2) * t83 / 0.2e1) * t83 + (t73 * mrSges(7,1) - t74 * mrSges(7,2) + Ifges(7,5) * t84 + Ifges(7,6) * t83 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t116 * mrSges(5,2) + t78 * mrSges(6,2) - t81 * mrSges(5,3) - t80 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t96) * t96 + (-t116 * mrSges(5,1) + t80 * mrSges(6,1) - t79 * mrSges(6,2) - t82 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t95 + (-Ifges(5,4) + Ifges(6,5)) * t96) * t95 + (t81 * mrSges(5,1) - t78 * mrSges(6,1) - t82 * mrSges(5,2) + t79 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t101 + (Ifges(6,4) + Ifges(5,5)) * t96 + (-Ifges(5,6) + Ifges(6,6)) * t95) * t101 + (t98 * (-mrSges(4,1) * t111 + mrSges(4,2) * t109) + (-t91 * t109 + t92 * t111) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t109 + Ifges(4,6) * t111) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t106 * mrSges(3,1) - t105 * mrSges(3,2) + (t105 ^ 2 + t106 ^ 2) * t121 * pkin(1)) * pkin(1) + Ifges(4,2) * t111 ^ 2 / 0.2e1 + (Ifges(4,4) * t111 + Ifges(4,1) * t109 / 0.2e1) * t109) * qJD(1)) * qJD(1);
T  = t1;
