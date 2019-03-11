% Calculate kinetic energy for
% S6RPPPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d5,d6,theta3,theta4]';
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
% Datum: 2019-03-09 01:34
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPPRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPPRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:33:14
% EndTime: 2019-03-09 01:33:15
% DurationCPUTime: 0.32s
% Computational Cost: add. (342->78), mult. (649->122), div. (0->0), fcn. (360->8), ass. (0->36)
t101 = sin(pkin(10));
t103 = cos(pkin(10));
t107 = sin(qJ(5));
t109 = cos(qJ(5));
t91 = (t101 * t107 - t103 * t109) * qJD(1);
t114 = m(3) / 0.2e1;
t100 = t103 * qJD(3);
t102 = sin(pkin(9));
t104 = cos(pkin(9));
t113 = qJ(2) * qJD(1);
t94 = qJD(2) + (-pkin(1) - pkin(2)) * qJD(1);
t89 = t102 * t94 + t104 * t113;
t87 = -qJD(1) * qJ(4) + t89;
t78 = t100 + (pkin(7) * qJD(1) - t87) * t101;
t112 = qJD(1) * t103;
t81 = t101 * qJD(3) + t103 * t87;
t79 = -pkin(7) * t112 + t81;
t74 = t107 * t78 + t109 * t79;
t88 = -t102 * t113 + t104 * t94;
t73 = -t107 * t79 + t109 * t78;
t86 = qJD(1) * pkin(3) + qJD(4) - t88;
t82 = pkin(4) * t112 + t86;
t108 = cos(qJ(6));
t106 = sin(qJ(6));
t98 = -qJD(1) * pkin(1) + qJD(2);
t92 = (-t101 * t109 - t103 * t107) * qJD(1);
t90 = qJD(6) - t91;
t85 = t106 * qJD(5) + t108 * t92;
t84 = t108 * qJD(5) - t106 * t92;
t80 = -t101 * t87 + t100;
t75 = -t91 * pkin(5) - t92 * pkin(8) + t82;
t72 = qJD(5) * pkin(8) + t74;
t71 = -qJD(5) * pkin(5) - t73;
t70 = t106 * t75 + t108 * t72;
t69 = -t106 * t72 + t108 * t75;
t1 = t98 ^ 2 * t114 + m(4) * (qJD(3) ^ 2 + t88 ^ 2 + t89 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t82 ^ 2) / 0.2e1 + m(5) * (t80 ^ 2 + t81 ^ 2 + t86 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + (t82 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,1) * t92 / 0.2e1) * t92 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t90 / 0.2e1) * t90 + (-t82 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t92 + Ifges(6,2) * t91 / 0.2e1) * t91 + (t71 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t90 + Ifges(7,1) * t85 / 0.2e1) * t85 + (-t71 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,6) * t90 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,5) * t92 + Ifges(6,6) * t91 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t86 * (mrSges(5,1) * t103 - mrSges(5,2) * t101) - t88 * mrSges(4,1) + t89 * mrSges(4,2) - t98 * mrSges(3,1) + (t80 * t101 - t81 * t103) * mrSges(5,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + (qJ(2) * t114 + mrSges(3,3)) * qJ(2) + t103 ^ 2 * Ifges(5,2) / 0.2e1 + (Ifges(5,4) * t103 + Ifges(5,1) * t101 / 0.2e1) * t101) * qJD(1)) * qJD(1);
T  = t1;
