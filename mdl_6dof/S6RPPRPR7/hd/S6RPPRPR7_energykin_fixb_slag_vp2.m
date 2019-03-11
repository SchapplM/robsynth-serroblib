% Calculate kinetic energy for
% S6RPPRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta3,theta5]';
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
% Datum: 2019-03-09 01:54
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR7_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:52:52
% EndTime: 2019-03-09 01:52:53
% DurationCPUTime: 0.41s
% Computational Cost: add. (486->87), mult. (1028->137), div. (0->0), fcn. (686->8), ass. (0->38)
t109 = cos(pkin(9));
t118 = t109 ^ 2;
t117 = m(3) / 0.2e1;
t106 = sin(pkin(10));
t108 = cos(pkin(10));
t111 = sin(qJ(4));
t113 = cos(qJ(4));
t107 = sin(pkin(9));
t99 = qJD(2) + (-pkin(1) - qJ(3)) * qJD(1);
t115 = -pkin(7) * qJD(1) + t99;
t92 = t115 * t107;
t93 = t115 * t109;
t86 = t111 * t93 + t113 * t92;
t83 = qJD(4) * qJ(5) + t86;
t95 = (t107 * t113 + t109 * t111) * qJD(1);
t96 = (-t107 * t111 + t109 * t113) * qJD(1);
t101 = qJD(1) * qJ(2) + qJD(3);
t97 = t107 * qJD(1) * pkin(3) + t101;
t84 = pkin(4) * t95 - qJ(5) * t96 + t97;
t75 = t106 * t84 + t108 * t83;
t116 = t107 ^ 2 + t118;
t74 = -t106 * t83 + t108 * t84;
t85 = -t111 * t92 + t113 * t93;
t82 = -qJD(4) * pkin(4) + qJD(5) - t85;
t112 = cos(qJ(6));
t110 = sin(qJ(6));
t102 = -qJD(1) * pkin(1) + qJD(2);
t94 = qJD(6) + t95;
t88 = qJD(4) * t106 + t108 * t96;
t87 = qJD(4) * t108 - t106 * t96;
t78 = t110 * t87 + t112 * t88;
t77 = -t110 * t88 + t112 * t87;
t76 = -pkin(5) * t87 + t82;
t73 = pkin(8) * t87 + t75;
t72 = pkin(5) * t95 - pkin(8) * t88 + t74;
t71 = t110 * t72 + t112 * t73;
t70 = -t110 * t73 + t112 * t72;
t1 = t102 ^ 2 * t117 + m(5) * (t85 ^ 2 + t86 ^ 2 + t97 ^ 2) / 0.2e1 + m(4) * (t116 * t99 ^ 2 + t101 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t82 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t76 ^ 2) / 0.2e1 + (t97 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t96 / 0.2e1) * t96 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t94 / 0.2e1) * t94 + (t82 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (-t82 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t88 + Ifges(6,2) * t87 / 0.2e1) * t87 + (t76 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t94 + Ifges(7,1) * t78 / 0.2e1) * t78 + (t85 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,5) * t96 + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (-t76 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t78 + Ifges(7,6) * t94 + Ifges(7,2) * t77 / 0.2e1) * t77 + (t97 * mrSges(5,1) + t74 * mrSges(6,1) - t75 * mrSges(6,2) - t86 * mrSges(5,3) - Ifges(5,4) * t96 + Ifges(6,5) * t88 - Ifges(5,6) * qJD(4) + Ifges(6,6) * t87 + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t95) * t95 + (t101 * (mrSges(4,1) * t107 + mrSges(4,2) * t109) + t102 * mrSges(3,2) - t116 * t99 * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + (qJ(2) * t117 + mrSges(3,3)) * qJ(2) + Ifges(4,1) * t118 / 0.2e1 + (-Ifges(4,4) * t109 + Ifges(4,2) * t107 / 0.2e1) * t107) * qJD(1)) * qJD(1);
T  = t1;
