% Calculate kinetic energy for
% S6PRPRRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d5,theta1,theta3]';
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
% Datum: 2019-03-08 20:12
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 20:09:09
% EndTime: 2019-03-08 20:09:10
% DurationCPUTime: 0.36s
% Computational Cost: add. (387->89), mult. (910->135), div. (0->0), fcn. (657->10), ass. (0->38)
t114 = sin(pkin(11));
t116 = cos(pkin(11));
t119 = sin(qJ(4));
t121 = cos(qJ(4));
t105 = (t114 * t119 - t116 * t121) * qJD(2);
t129 = cos(qJ(5));
t118 = sin(qJ(5));
t120 = sin(qJ(2));
t115 = sin(pkin(6));
t127 = qJD(1) * t115;
t109 = qJD(2) * qJ(3) + t120 * t127;
t117 = cos(pkin(6));
t126 = qJD(1) * t117;
t111 = t116 * t126;
t128 = pkin(8) * qJD(2);
t97 = t111 + (-t109 - t128) * t114;
t102 = t116 * t109 + t114 * t126;
t98 = t116 * t128 + t102;
t91 = t119 * t97 + t121 * t98;
t89 = qJD(4) * pkin(9) + t91;
t122 = cos(qJ(2));
t125 = -t122 * t127 + qJD(3);
t103 = (-pkin(3) * t116 - pkin(2)) * qJD(2) + t125;
t106 = (t114 * t121 + t116 * t119) * qJD(2);
t93 = t105 * pkin(4) - t106 * pkin(9) + t103;
t85 = t118 * t93 + t129 * t89;
t90 = -t119 * t98 + t121 * t97;
t88 = -qJD(4) * pkin(4) - t90;
t84 = -t118 * t89 + t129 * t93;
t108 = -qJD(2) * pkin(2) + t125;
t104 = qJD(5) + t105;
t101 = -t114 * t109 + t111;
t100 = t118 * qJD(4) + t106 * t129;
t99 = -qJD(4) * t129 + t118 * t106;
t86 = t99 * pkin(5) - t100 * qJ(6) + t88;
t83 = t104 * qJ(6) + t85;
t82 = -t104 * pkin(5) + qJD(6) - t84;
t1 = m(7) * (t82 ^ 2 + t83 ^ 2 + t86 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t88 ^ 2) / 0.2e1 + m(4) * (t101 ^ 2 + t102 ^ 2 + t108 ^ 2) / 0.2e1 + m(5) * (t103 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + (t103 * mrSges(5,2) - t90 * mrSges(5,3) + Ifges(5,1) * t106 / 0.2e1) * t106 - (-t103 * mrSges(5,1) + t91 * mrSges(5,3) + Ifges(5,4) * t106 - Ifges(5,2) * t105 / 0.2e1) * t105 + (m(2) / 0.2e1 + m(3) * (t117 ^ 2 + (t120 ^ 2 + t122 ^ 2) * t115 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t90 * mrSges(5,1) - t91 * mrSges(5,2) + Ifges(5,5) * t106 - Ifges(5,6) * t105 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t88 * mrSges(6,1) + t86 * mrSges(7,1) - t83 * mrSges(7,2) - t85 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t99) * t99 + (t84 * mrSges(6,1) - t82 * mrSges(7,1) - t85 * mrSges(6,2) + t83 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t104 + (-Ifges(6,6) + Ifges(7,6)) * t99) * t104 + (t88 * mrSges(6,2) + t82 * mrSges(7,2) - t84 * mrSges(6,3) - t86 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t100 + (-Ifges(6,4) + Ifges(7,5)) * t99 + (Ifges(7,4) + Ifges(6,5)) * t104) * t100 + (t108 * (-mrSges(4,1) * t116 + mrSges(4,2) * t114) + (Ifges(3,3) / 0.2e1 + Ifges(4,2) * t116 ^ 2 / 0.2e1 + (Ifges(4,4) * t116 + Ifges(4,1) * t114 / 0.2e1) * t114) * qJD(2) + (mrSges(3,1) * t122 - mrSges(3,2) * t120) * t127 + (-t101 * t114 + t102 * t116) * mrSges(4,3)) * qJD(2);
T  = t1;
