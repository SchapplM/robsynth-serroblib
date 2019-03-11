% Calculate kinetic energy for
% S6PRPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 19:23
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:21:21
% EndTime: 2019-03-08 19:21:22
% DurationCPUTime: 0.22s
% Computational Cost: add. (233->69), mult. (436->111), div. (0->0), fcn. (239->10), ass. (0->35)
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t115 = cos(qJ(2));
t106 = sin(pkin(6));
t121 = qJD(1) * t106;
t118 = -t115 * t121 + qJD(3);
t93 = (-pkin(2) - pkin(3)) * qJD(2) + t118;
t112 = sin(qJ(2));
t99 = qJD(2) * qJ(3) + t112 * t121;
t91 = t105 * t93 + t107 * t99;
t108 = cos(pkin(6));
t101 = -t108 * qJD(1) + qJD(4);
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t89 = -qJD(2) * pkin(8) + t91;
t85 = t111 * t101 + t114 * t89;
t120 = qJD(2) * t111;
t119 = t114 * qJD(2);
t90 = -t105 * t99 + t107 * t93;
t88 = qJD(2) * pkin(4) - t90;
t84 = t114 * t101 - t111 * t89;
t117 = qJD(1) ^ 2;
t113 = cos(qJ(6));
t110 = sin(qJ(6));
t103 = t108 ^ 2 * t117;
t102 = qJD(6) + t119;
t98 = t110 * qJD(5) - t113 * t120;
t97 = t113 * qJD(5) + t110 * t120;
t96 = -qJD(2) * pkin(2) + t118;
t86 = (pkin(5) * t114 + pkin(9) * t111) * qJD(2) + t88;
t83 = qJD(5) * pkin(9) + t85;
t82 = -qJD(5) * pkin(5) - t84;
t81 = t110 * t86 + t113 * t83;
t80 = -t110 * t83 + t113 * t86;
t1 = m(3) * (t103 + (t112 ^ 2 + t115 ^ 2) * t117 * t106 ^ 2) / 0.2e1 + m(2) * t117 / 0.2e1 + m(5) * (t101 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(4) * (t96 ^ 2 + t99 ^ 2 + t103) / 0.2e1 + m(7) * (t80 ^ 2 + t81 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t84 ^ 2 + t85 ^ 2 + t88 ^ 2) / 0.2e1 + (t82 * mrSges(7,2) - t80 * mrSges(7,3) + Ifges(7,1) * t98 / 0.2e1) * t98 + (t84 * mrSges(6,1) - t85 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t82 * mrSges(7,1) + t81 * mrSges(7,3) + Ifges(7,4) * t98 + Ifges(7,2) * t97 / 0.2e1) * t97 + (t80 * mrSges(7,1) - t81 * mrSges(7,2) + Ifges(7,5) * t98 + Ifges(7,6) * t97 + Ifges(7,3) * t102 / 0.2e1) * t102 + (-t96 * mrSges(4,1) - t90 * mrSges(5,1) + t91 * mrSges(5,2) + t99 * mrSges(4,3) + (mrSges(3,1) * t115 - mrSges(3,2) * t112) * t121 + (t88 * mrSges(6,1) - t85 * mrSges(6,3) - Ifges(6,6) * qJD(5) + Ifges(6,2) * t119 / 0.2e1) * t114 + (-t88 * mrSges(6,2) + t84 * mrSges(6,3) - Ifges(6,5) * qJD(5)) * t111 + (Ifges(3,3) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + (Ifges(6,4) * t114 + Ifges(6,1) * t111 / 0.2e1) * t111) * qJD(2)) * qJD(2);
T  = t1;
