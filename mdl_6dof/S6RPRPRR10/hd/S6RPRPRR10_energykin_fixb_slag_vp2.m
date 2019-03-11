% Calculate kinetic energy for
% S6RPRPRR10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta4]';
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
% Datum: 2019-03-09 04:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR10_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 04:07:33
% EndTime: 2019-03-09 04:07:33
% DurationCPUTime: 0.43s
% Computational Cost: add. (561->95), mult. (1128->149), div. (0->0), fcn. (740->8), ass. (0->39)
t108 = sin(qJ(5));
t111 = cos(qJ(5));
t109 = sin(qJ(3));
t103 = t109 * qJD(1);
t105 = sin(pkin(10));
t106 = cos(pkin(10));
t112 = cos(qJ(3));
t94 = (pkin(3) * t109 - qJ(4) * t112 + qJ(2)) * qJD(1);
t100 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t95 = qJD(3) * qJ(4) + t109 * t100;
t85 = -t105 * t95 + t106 * t94;
t114 = qJD(1) * t112;
t97 = qJD(3) * t105 + t106 * t114;
t82 = pkin(4) * t103 - pkin(8) * t97 + t85;
t86 = t105 * t94 + t106 * t95;
t96 = qJD(3) * t106 - t105 * t114;
t84 = pkin(8) * t96 + t86;
t76 = t108 * t82 + t111 * t84;
t115 = t100 * mrSges(4,3);
t101 = t103 + qJD(5);
t75 = -t108 * t84 + t111 * t82;
t93 = -qJD(3) * pkin(3) - t112 * t100 + qJD(4);
t89 = -pkin(4) * t96 + t93;
t113 = qJD(1) ^ 2;
t110 = cos(qJ(6));
t107 = sin(qJ(6));
t104 = t113 * qJ(2) ^ 2;
t102 = -qJD(1) * pkin(1) + qJD(2);
t99 = qJD(6) + t101;
t88 = t108 * t96 + t111 * t97;
t87 = -t108 * t97 + t111 * t96;
t79 = -pkin(5) * t87 + t89;
t78 = t107 * t87 + t110 * t88;
t77 = -t107 * t88 + t110 * t87;
t74 = pkin(9) * t87 + t76;
t73 = pkin(5) * t101 - pkin(9) * t88 + t75;
t72 = t107 * t73 + t110 * t74;
t71 = -t107 * t74 + t110 * t73;
t1 = m(3) * (t102 ^ 2 + t104) / 0.2e1 + m(4) * (t104 + (t109 ^ 2 + t112 ^ 2) * t100 ^ 2) / 0.2e1 + m(7) * (t71 ^ 2 + t72 ^ 2 + t79 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t93 ^ 2) / 0.2e1 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,3) * t99 / 0.2e1) * t99 + (t93 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,1) * t97 / 0.2e1) * t97 + (t89 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,1) * t88 / 0.2e1) * t88 + (qJ(2) * mrSges(3,3) + Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1) * t113 + (Ifges(4,3) * qJD(3) / 0.2e1 + (t112 * mrSges(4,1) - t109 * mrSges(4,2)) * t100) * qJD(3) + (-t89 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t88 + Ifges(6,2) * t87 / 0.2e1) * t87 + (t79 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,5) * t99 + Ifges(7,1) * t78 / 0.2e1) * t78 + (t102 * mrSges(3,2) + (qJ(2) * mrSges(4,2) * qJD(1) + Ifges(4,5) * qJD(3) + (-t115 + Ifges(4,1) * qJD(1) / 0.2e1) * t112) * t112 + (Ifges(5,5) * t97 - Ifges(4,6) * qJD(3) - t86 * mrSges(5,2) + t85 * mrSges(5,1) + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t112) * qJD(1) + (-t115 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(1)) * t109) * t109) * qJD(1) + (Ifges(5,6) * t103 - t93 * mrSges(5,1) + t86 * mrSges(5,3) + Ifges(5,4) * t97 + Ifges(5,2) * t96 / 0.2e1) * t96 + (-t79 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t78 + Ifges(7,6) * t99 + Ifges(7,2) * t77 / 0.2e1) * t77 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,5) * t88 + Ifges(6,6) * t87 + Ifges(6,3) * t101 / 0.2e1) * t101;
T  = t1;
