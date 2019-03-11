% Calculate kinetic energy for
% S6PRPRRP1
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
% Datum: 2019-03-08 19:58
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPRRP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:55:08
% EndTime: 2019-03-08 19:55:08
% DurationCPUTime: 0.30s
% Computational Cost: add. (275->80), mult. (597->119), div. (0->0), fcn. (386->10), ass. (0->34)
t111 = sin(qJ(5));
t114 = cos(qJ(5));
t110 = cos(pkin(6));
t104 = qJD(1) * t110 + qJD(3);
t112 = sin(qJ(4));
t115 = cos(qJ(4));
t116 = cos(qJ(2));
t108 = sin(pkin(6));
t122 = qJD(1) * t108;
t100 = qJD(2) * pkin(2) + t116 * t122;
t107 = sin(pkin(11));
t109 = cos(pkin(11));
t113 = sin(qJ(2));
t119 = t113 * t122;
t96 = t107 * t100 + t109 * t119;
t94 = qJD(2) * pkin(8) + t96;
t87 = t112 * t104 + t115 * t94;
t85 = qJD(4) * pkin(9) + t87;
t95 = t100 * t109 - t107 * t119;
t90 = (-pkin(4) * t115 - pkin(9) * t112 - pkin(3)) * qJD(2) - t95;
t81 = t111 * t90 + t114 * t85;
t121 = t112 * qJD(2);
t120 = t115 * qJD(2);
t80 = -t111 * t85 + t114 * t90;
t86 = t104 * t115 - t112 * t94;
t84 = -qJD(4) * pkin(4) - t86;
t105 = qJD(5) - t120;
t99 = qJD(4) * t111 + t114 * t121;
t98 = qJD(4) * t114 - t111 * t121;
t93 = -qJD(2) * pkin(3) - t95;
t82 = -pkin(5) * t98 + qJD(6) + t84;
t79 = qJ(6) * t98 + t81;
t78 = pkin(5) * t105 - qJ(6) * t99 + t80;
t1 = m(5) * (t86 ^ 2 + t87 ^ 2 + t93 ^ 2) / 0.2e1 + m(7) * (t78 ^ 2 + t79 ^ 2 + t82 ^ 2) / 0.2e1 + m(6) * (t80 ^ 2 + t81 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t104 ^ 2 + t95 ^ 2 + t96 ^ 2) / 0.2e1 + (t86 * mrSges(5,1) - t87 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (m(2) / 0.2e1 + m(3) * (t110 ^ 2 + (t113 ^ 2 + t116 ^ 2) * t108 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t84 * mrSges(6,2) + t82 * mrSges(7,2) - t80 * mrSges(6,3) - t78 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t99) * t99 + (-t84 * mrSges(6,1) - t82 * mrSges(7,1) + t81 * mrSges(6,3) + t79 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t98 + (Ifges(6,4) + Ifges(7,4)) * t99) * t98 + (t80 * mrSges(6,1) + t78 * mrSges(7,1) - t81 * mrSges(6,2) - t79 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t105 + (Ifges(6,5) + Ifges(7,5)) * t99 + (Ifges(6,6) + Ifges(7,6)) * t98) * t105 + (t95 * mrSges(4,1) - t96 * mrSges(4,2) + (mrSges(3,1) * t116 - mrSges(3,2) * t113) * t122 + (-t93 * mrSges(5,1) + t87 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t120 / 0.2e1) * t115 + (t93 * mrSges(5,2) - t86 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t112 + (Ifges(4,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (Ifges(5,4) * t115 + Ifges(5,1) * t112 / 0.2e1) * t112) * qJD(2)) * qJD(2);
T  = t1;
