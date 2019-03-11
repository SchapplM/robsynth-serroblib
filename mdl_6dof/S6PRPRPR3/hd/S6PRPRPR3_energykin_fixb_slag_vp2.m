% Calculate kinetic energy for
% S6PRPRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d4,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:38
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:34:21
% EndTime: 2019-03-08 19:34:21
% DurationCPUTime: 0.29s
% Computational Cost: add. (233->82), mult. (507->120), div. (0->0), fcn. (304->10), ass. (0->37)
t123 = -pkin(4) - pkin(9);
t105 = sin(pkin(11));
t107 = cos(pkin(11));
t111 = sin(qJ(2));
t106 = sin(pkin(6));
t122 = qJD(1) * t106;
t119 = t111 * t122;
t114 = cos(qJ(2));
t98 = qJD(2) * pkin(2) + t114 * t122;
t94 = t105 * t98 + t107 * t119;
t108 = cos(pkin(6));
t102 = qJD(1) * t108 + qJD(3);
t110 = sin(qJ(4));
t113 = cos(qJ(4));
t92 = qJD(2) * pkin(8) + t94;
t87 = t110 * t102 + t113 * t92;
t121 = qJD(2) * t110;
t120 = qJD(2) * t113;
t118 = -qJ(5) * t110 - pkin(3);
t93 = -t105 * t119 + t107 * t98;
t86 = t102 * t113 - t110 * t92;
t84 = -qJD(4) * qJ(5) - t87;
t117 = qJD(5) - t86;
t112 = cos(qJ(6));
t109 = sin(qJ(6));
t103 = qJD(6) + t121;
t97 = qJD(4) * t112 - t109 * t120;
t96 = -qJD(4) * t109 - t112 * t120;
t91 = -qJD(2) * pkin(3) - t93;
t88 = (-pkin(4) * t113 + t118) * qJD(2) - t93;
t85 = (t123 * t113 + t118) * qJD(2) - t93;
t83 = -qJD(4) * pkin(4) + t117;
t82 = pkin(5) * t120 - t84;
t81 = pkin(5) * t121 + t123 * qJD(4) + t117;
t80 = t109 * t81 + t112 * t85;
t79 = -t109 * t85 + t112 * t81;
t1 = m(4) * (t102 ^ 2 + t93 ^ 2 + t94 ^ 2) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t83 ^ 2 + t84 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t91 ^ 2) / 0.2e1 + (t82 * mrSges(7,2) - t79 * mrSges(7,3) + Ifges(7,1) * t97 / 0.2e1) * t97 + (-t82 * mrSges(7,1) + t80 * mrSges(7,3) + Ifges(7,4) * t97 + Ifges(7,2) * t96 / 0.2e1) * t96 + (m(2) / 0.2e1 + m(3) * (t108 ^ 2 + (t111 ^ 2 + t114 ^ 2) * t106 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t79 * mrSges(7,1) - t80 * mrSges(7,2) + Ifges(7,5) * t97 + Ifges(7,6) * t96 + Ifges(7,3) * t103 / 0.2e1) * t103 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + t83 * mrSges(6,2) - t84 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4)) * qJD(4) + (t93 * mrSges(4,1) - t94 * mrSges(4,2) + (mrSges(3,1) * t114 - mrSges(3,2) * t111) * t122 + (-t91 * mrSges(5,1) - t84 * mrSges(6,1) + t88 * mrSges(6,2) + t87 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t120 + (-Ifges(6,5) + Ifges(5,6)) * qJD(4)) * t113 + (t83 * mrSges(6,1) + t91 * mrSges(5,2) - t86 * mrSges(5,3) - t88 * mrSges(6,3) + (-Ifges(6,4) + Ifges(5,5)) * qJD(4)) * t110 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + ((Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t110 + (Ifges(5,4) + Ifges(6,6)) * t113) * t110) * qJD(2)) * qJD(2);
T  = t1;
