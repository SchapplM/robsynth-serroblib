% Calculate kinetic energy for
% S6RPPRRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d5,d6,theta2]';
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
% Datum: 2019-03-09 02:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRRR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:22:54
% EndTime: 2019-03-09 02:22:55
% DurationCPUTime: 0.37s
% Computational Cost: add. (316->82), mult. (600->126), div. (0->0), fcn. (328->8), ass. (0->36)
t100 = sin(pkin(10));
t111 = -pkin(1) * t100 - qJ(3);
t93 = t111 * qJD(1);
t115 = t93 ^ 2;
t114 = m(3) / 0.2e1;
t103 = sin(qJ(5));
t106 = cos(qJ(5));
t104 = sin(qJ(4));
t107 = cos(qJ(4));
t101 = cos(pkin(10));
t112 = -pkin(1) * t101 - pkin(2);
t89 = qJD(3) + (-pkin(7) + t112) * qJD(1);
t86 = t107 * qJD(2) + t104 * t89;
t82 = qJD(4) * pkin(8) + t86;
t87 = (pkin(4) * t104 - pkin(8) * t107 - t111) * qJD(1);
t76 = t103 * t87 + t106 * t82;
t96 = t104 * qJD(1) + qJD(5);
t113 = qJD(1) * t107;
t75 = -t103 * t82 + t106 * t87;
t85 = -t104 * qJD(2) + t107 * t89;
t81 = -qJD(4) * pkin(4) - t85;
t108 = qJD(2) ^ 2;
t105 = cos(qJ(6));
t102 = sin(qJ(6));
t95 = qJD(6) + t96;
t92 = t112 * qJD(1) + qJD(3);
t91 = qJD(4) * t103 + t106 * t113;
t90 = qJD(4) * t106 - t103 * t113;
t79 = t102 * t90 + t105 * t91;
t78 = -t102 * t91 + t105 * t90;
t77 = -pkin(5) * t90 + t81;
t74 = pkin(9) * t90 + t76;
t73 = pkin(5) * t96 - pkin(9) * t91 + t75;
t72 = t102 * t73 + t105 * t74;
t71 = -t102 * t74 + t105 * t73;
t1 = m(4) * (t92 ^ 2 + t108 + t115) / 0.2e1 + t108 * t114 + m(5) * (t85 ^ 2 + t86 ^ 2 + t115) / 0.2e1 + m(7) * (t71 ^ 2 + t72 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t81 ^ 2) / 0.2e1 + (t75 * mrSges(6,1) - t76 * mrSges(6,2) + Ifges(6,3) * t96 / 0.2e1) * t96 + (t71 * mrSges(7,1) - t72 * mrSges(7,2) + Ifges(7,3) * t95 / 0.2e1) * t95 + (t85 * mrSges(5,1) - t86 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t81 * mrSges(6,2) - t75 * mrSges(6,3) + Ifges(6,5) * t96 + Ifges(6,1) * t91 / 0.2e1) * t91 + (t77 * mrSges(7,2) - t71 * mrSges(7,3) + Ifges(7,5) * t95 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t81 * mrSges(6,1) + t76 * mrSges(6,3) + Ifges(6,4) * t91 + Ifges(6,6) * t96 + Ifges(6,2) * t90 / 0.2e1) * t90 + (-t77 * mrSges(7,1) + t72 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t95 + Ifges(7,2) * t78 / 0.2e1) * t78 + (t92 * mrSges(4,2) - t93 * mrSges(4,3) + (-t93 * mrSges(5,2) - t85 * mrSges(5,3) + Ifges(5,5) * qJD(4) + Ifges(5,1) * t113 / 0.2e1) * t107 + (-t93 * mrSges(5,1) - t86 * mrSges(5,3) - Ifges(5,6) * qJD(4)) * t104 + (Ifges(3,3) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(2,3) / 0.2e1 + (mrSges(3,1) * t101 - mrSges(3,2) * t100 + (t100 ^ 2 + t101 ^ 2) * t114 * pkin(1)) * pkin(1) + (-Ifges(5,4) * t107 + Ifges(5,2) * t104 / 0.2e1) * t104) * qJD(1)) * qJD(1);
T  = t1;
