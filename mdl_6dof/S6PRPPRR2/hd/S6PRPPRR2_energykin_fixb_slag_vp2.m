% Calculate kinetic energy for
% S6PRPPRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d5,d6,theta1,theta3]';
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
% Datum: 2019-03-08 19:20
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRPPRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(11,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRPPRR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 19:17:36
% EndTime: 2019-03-08 19:17:36
% DurationCPUTime: 0.25s
% Computational Cost: add. (202->69), mult. (425->108), div. (0->0), fcn. (255->10), ass. (0->34)
t102 = sin(pkin(11));
t104 = cos(pkin(11));
t108 = sin(qJ(2));
t103 = sin(pkin(6));
t117 = qJD(1) * t103;
t115 = t108 * t117;
t111 = cos(qJ(2));
t95 = qJD(2) * pkin(2) + t111 * t117;
t92 = t102 * t95 + t104 * t115;
t89 = qJD(2) * qJ(4) + t92;
t118 = t89 ^ 2;
t107 = sin(qJ(5));
t110 = cos(qJ(5));
t91 = -t102 * t115 + t104 * t95;
t114 = qJD(4) - t91;
t87 = (-pkin(3) - pkin(8)) * qJD(2) + t114;
t105 = cos(pkin(6));
t99 = qJD(1) * t105 + qJD(3);
t84 = t107 * t87 + t110 * t99;
t116 = qJD(2) * t110;
t83 = -t107 * t99 + t110 * t87;
t109 = cos(qJ(6));
t106 = sin(qJ(6));
t100 = qJD(2) * t107 + qJD(6);
t98 = t99 ^ 2;
t94 = qJD(5) * t106 + t109 * t116;
t93 = qJD(5) * t109 - t106 * t116;
t88 = -qJD(2) * pkin(3) + t114;
t85 = (pkin(5) * t107 - pkin(9) * t110 + qJ(4)) * qJD(2) + t92;
t82 = qJD(5) * pkin(9) + t84;
t81 = -qJD(5) * pkin(5) - t83;
t80 = t106 * t85 + t109 * t82;
t79 = -t106 * t82 + t109 * t85;
t1 = m(6) * (t83 ^ 2 + t84 ^ 2 + t118) / 0.2e1 + m(5) * (t88 ^ 2 + t118 + t98) / 0.2e1 + m(4) * (t91 ^ 2 + t92 ^ 2 + t98) / 0.2e1 + m(7) * (t79 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + (t81 * mrSges(7,2) - t79 * mrSges(7,3) + Ifges(7,1) * t94 / 0.2e1) * t94 + (t83 * mrSges(6,1) - t84 * mrSges(6,2) + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t81 * mrSges(7,1) + t80 * mrSges(7,3) + Ifges(7,4) * t94 + Ifges(7,2) * t93 / 0.2e1) * t93 + (m(2) / 0.2e1 + m(3) * (t105 ^ 2 + (t108 ^ 2 + t111 ^ 2) * t103 ^ 2) / 0.2e1) * qJD(1) ^ 2 + (t79 * mrSges(7,1) - t80 * mrSges(7,2) + Ifges(7,5) * t94 + Ifges(7,6) * t93 + Ifges(7,3) * t100 / 0.2e1) * t100 + (t91 * mrSges(4,1) - t92 * mrSges(4,2) + t88 * mrSges(5,2) + t89 * mrSges(5,3) + (mrSges(3,1) * t111 - mrSges(3,2) * t108) * t117 + (t89 * mrSges(6,2) - t83 * mrSges(6,3) + Ifges(6,5) * qJD(5) + Ifges(6,1) * t116 / 0.2e1) * t110 + (t89 * mrSges(6,1) - t84 * mrSges(6,3) - Ifges(6,6) * qJD(5)) * t107 + (Ifges(5,1) / 0.2e1 + Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (-Ifges(6,4) * t110 + Ifges(6,2) * t107 / 0.2e1) * t107) * qJD(2)) * qJD(2);
T  = t1;
