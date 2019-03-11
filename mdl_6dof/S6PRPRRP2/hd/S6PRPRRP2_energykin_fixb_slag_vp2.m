% Calculate kinetic energy for
% S6PRPRRP2
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
% Datum: 2019-03-08 20:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:59:51
% EndTime: 2019-03-08 19:59:51
% DurationCPUTime: 0.30s
% Computational Cost: add. (275->80), mult. (593->119), div. (0->0), fcn. (382->10), ass. (0->34)
t122 = cos(qJ(5));
t111 = sin(qJ(5));
t110 = cos(pkin(6));
t103 = qJD(1) * t110 + qJD(3);
t112 = sin(qJ(4));
t114 = cos(qJ(4));
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t113 = sin(qJ(2));
t108 = sin(pkin(6));
t121 = qJD(1) * t108;
t118 = t113 * t121;
t115 = cos(qJ(2));
t99 = qJD(2) * pkin(2) + t115 * t121;
t95 = t107 * t99 + t109 * t118;
t93 = qJD(2) * pkin(8) + t95;
t87 = t112 * t103 + t114 * t93;
t85 = qJD(4) * pkin(9) + t87;
t94 = -t107 * t118 + t109 * t99;
t89 = (-pkin(4) * t114 - pkin(9) * t112 - pkin(3)) * qJD(2) - t94;
t81 = t111 * t89 + t122 * t85;
t120 = qJD(2) * t114;
t119 = t112 * qJD(2);
t86 = t103 * t114 - t112 * t93;
t84 = -qJD(4) * pkin(4) - t86;
t80 = -t111 * t85 + t122 * t89;
t104 = qJD(5) - t120;
t98 = t111 * qJD(4) + t119 * t122;
t97 = -qJD(4) * t122 + t111 * t119;
t92 = -qJD(2) * pkin(3) - t94;
t82 = pkin(5) * t97 - qJ(6) * t98 + t84;
t79 = t104 * qJ(6) + t81;
t78 = -t104 * pkin(5) + qJD(6) - t80;
t1 = m(4) * (t103 ^ 2 + t94 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t92 ^ 2) / 0.2e1 + m(7) * (t78 ^ 2 + t79 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (m(2) / 0.2e1 + m(3) * (t110 ^ 2 + (t113 ^ 2 + t115 ^ 2) * t108 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t84 * mrSges(6,2) + t78 * mrSges(7,2) - t80 * mrSges(6,3) - t82 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t98) * t98 + (t84 * mrSges(6,1) + t82 * mrSges(7,1) - t79 * mrSges(7,2) - t81 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t97 + (-Ifges(6,4) + Ifges(7,5)) * t98) * t97 + (t80 * mrSges(6,1) - t78 * mrSges(7,1) - t81 * mrSges(6,2) + t79 * mrSges(7,3) + (Ifges(6,3) / 0.2e1 + Ifges(7,2) / 0.2e1) * t104 + (Ifges(7,4) + Ifges(6,5)) * t98 + (-Ifges(6,6) + Ifges(7,6)) * t97) * t104 + (t94 * mrSges(4,1) - t95 * mrSges(4,2) + (mrSges(3,1) * t115 - mrSges(3,2) * t113) * t121 + (-t92 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t120 / 0.2e1) * t114 + (t92 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t112 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t114 + Ifges(5,1) * t112 / 0.2e1) * t112) * qJD(2)) * qJD(2);
T  = t1;
