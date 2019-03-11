% Calculate kinetic energy for
% S6PRRPRP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d5,theta1]';
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
% Datum: 2019-03-08 21:44
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6PRRPRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRP4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-08 21:40:06
% EndTime: 2019-03-08 21:40:06
% DurationCPUTime: 0.29s
% Computational Cost: add. (266->92), mult. (557->125), div. (0->0), fcn. (317->8), ass. (0->35)
t118 = -pkin(3) - pkin(9);
t104 = sin(qJ(5));
t107 = cos(qJ(5));
t105 = sin(qJ(3));
t108 = cos(qJ(3));
t103 = cos(pkin(6));
t116 = qJD(1) * t103;
t106 = sin(qJ(2));
t102 = sin(pkin(6));
t117 = qJD(1) * t102;
t95 = qJD(2) * pkin(8) + t106 * t117;
t88 = -t105 * t95 + t108 * t116;
t111 = qJD(4) - t88;
t114 = t105 * qJD(2);
t82 = pkin(4) * t114 + t118 * qJD(3) + t111;
t112 = -qJ(4) * t105 - pkin(2);
t109 = cos(qJ(2));
t113 = t109 * t117;
t85 = -t113 + (t118 * t108 + t112) * qJD(2);
t78 = t104 * t82 + t107 * t85;
t89 = t105 * t116 + t108 * t95;
t115 = qJD(2) * t108;
t87 = -qJD(3) * qJ(4) - t89;
t77 = -t104 * t85 + t107 * t82;
t83 = pkin(4) * t115 - t87;
t98 = qJD(5) + t114;
t96 = -qJD(2) * pkin(2) - t113;
t94 = t107 * qJD(3) - t104 * t115;
t93 = -t104 * qJD(3) - t107 * t115;
t90 = -t113 + (-pkin(3) * t108 + t112) * qJD(2);
t86 = -qJD(3) * pkin(3) + t111;
t79 = -t93 * pkin(5) + qJD(6) + t83;
t76 = t93 * qJ(6) + t78;
t75 = t98 * pkin(5) - t94 * qJ(6) + t77;
t1 = m(4) * (t88 ^ 2 + t89 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t86 ^ 2 + t87 ^ 2 + t90 ^ 2) / 0.2e1 + m(7) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t77 ^ 2 + t78 ^ 2 + t83 ^ 2) / 0.2e1 + (m(3) * (t103 ^ 2 + (t106 ^ 2 + t109 ^ 2) * t102 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t77 * mrSges(6,1) + t75 * mrSges(7,1) - t78 * mrSges(6,2) - t76 * mrSges(7,2) + (Ifges(6,3) / 0.2e1 + Ifges(7,3) / 0.2e1) * t98) * t98 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + t86 * mrSges(5,2) - t87 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3)) * qJD(3) + (t83 * mrSges(6,2) + t79 * mrSges(7,2) - t77 * mrSges(6,3) - t75 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t94 + (Ifges(6,5) + Ifges(7,5)) * t98) * t94 + (-t83 * mrSges(6,1) - t79 * mrSges(7,1) + t78 * mrSges(6,3) + t76 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t93 + (Ifges(6,6) + Ifges(7,6)) * t98 + (Ifges(6,4) + Ifges(7,4)) * t94) * t93 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t109 - mrSges(3,2) * t106) * t117 + (-t96 * mrSges(4,1) - t87 * mrSges(5,1) + t90 * mrSges(5,2) + t89 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t115 + (-Ifges(5,5) + Ifges(4,6)) * qJD(3)) * t108 + (t86 * mrSges(5,1) + t96 * mrSges(4,2) - t88 * mrSges(4,3) - t90 * mrSges(5,3) + (-Ifges(5,4) + Ifges(4,5)) * qJD(3) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t105 + (Ifges(4,4) + Ifges(5,6)) * t108) * qJD(2)) * t105) * qJD(2);
T  = t1;
