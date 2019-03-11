% Calculate kinetic energy for
% S6RPRPRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-03-09 03:03
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:01:16
% EndTime: 2019-03-09 03:01:16
% DurationCPUTime: 0.44s
% Computational Cost: add. (408->94), mult. (903->138), div. (0->0), fcn. (574->8), ass. (0->35)
t100 = cos(pkin(10));
t103 = sin(qJ(3));
t105 = cos(qJ(3));
t98 = sin(pkin(10));
t91 = (t100 * t105 - t103 * t98) * qJD(1);
t111 = m(3) / 0.2e1;
t102 = sin(qJ(5));
t104 = cos(qJ(5));
t110 = qJ(4) * qJD(1);
t99 = sin(pkin(9));
t94 = (pkin(1) * t99 + pkin(7)) * qJD(1);
t97 = t105 * qJD(2);
t83 = qJD(3) * pkin(3) + t97 + (-t94 - t110) * t103;
t88 = t103 * qJD(2) + t105 * t94;
t86 = t105 * t110 + t88;
t76 = t100 * t86 + t98 * t83;
t74 = qJD(3) * pkin(8) + t76;
t101 = cos(pkin(9));
t109 = -pkin(1) * t101 - pkin(2);
t90 = qJD(4) + (-pkin(3) * t105 + t109) * qJD(1);
t92 = (t100 * t103 + t105 * t98) * qJD(1);
t79 = -pkin(4) * t91 - pkin(8) * t92 + t90;
t70 = t102 * t79 + t104 * t74;
t75 = t100 * t83 - t98 * t86;
t69 = -t102 * t74 + t104 * t79;
t73 = -qJD(3) * pkin(4) - t75;
t95 = t109 * qJD(1);
t89 = qJD(5) - t91;
t87 = -t103 * t94 + t97;
t85 = qJD(3) * t102 + t104 * t92;
t84 = qJD(3) * t104 - t102 * t92;
t71 = -pkin(5) * t84 + qJD(6) + t73;
t68 = qJ(6) * t84 + t70;
t67 = pkin(5) * t89 - qJ(6) * t85 + t69;
t1 = qJD(2) ^ 2 * t111 + m(5) * (t75 ^ 2 + t76 ^ 2 + t90 ^ 2) / 0.2e1 + m(4) * (t87 ^ 2 + t88 ^ 2 + t95 ^ 2) / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + (t90 * mrSges(5,2) - t75 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 + (-t90 * mrSges(5,1) + t76 * mrSges(5,3) + Ifges(5,4) * t92 + Ifges(5,2) * t91 / 0.2e1) * t91 + (t69 * mrSges(6,1) + t67 * mrSges(7,1) - t70 * mrSges(6,2) - t68 * mrSges(7,2) + (Ifges(7,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t89) * t89 + (t73 * mrSges(6,2) + t71 * mrSges(7,2) - t69 * mrSges(6,3) - t67 * mrSges(7,3) + (Ifges(6,1) / 0.2e1 + Ifges(7,1) / 0.2e1) * t85 + (Ifges(6,5) + Ifges(7,5)) * t89) * t85 + (-t73 * mrSges(6,1) - t71 * mrSges(7,1) + t70 * mrSges(6,3) + t68 * mrSges(7,3) + (Ifges(6,2) / 0.2e1 + Ifges(7,2) / 0.2e1) * t84 + (Ifges(6,6) + Ifges(7,6)) * t89 + (Ifges(6,4) + Ifges(7,4)) * t85) * t84 + (t87 * mrSges(4,1) + t75 * mrSges(5,1) - t88 * mrSges(4,2) - t76 * mrSges(5,2) + Ifges(5,5) * t92 + Ifges(5,6) * t91 + (Ifges(5,3) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3)) * qJD(3) + (t95 * (-mrSges(4,1) * t105 + mrSges(4,2) * t103) + (-t87 * t103 + t88 * t105) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t103 + Ifges(4,6) * t105) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t101 * mrSges(3,1) - t99 * mrSges(3,2) + (t101 ^ 2 + t99 ^ 2) * t111 * pkin(1)) * pkin(1) + Ifges(4,2) * t105 ^ 2 / 0.2e1 + (Ifges(4,4) * t105 + Ifges(4,1) * t103 / 0.2e1) * t103) * qJD(1)) * qJD(1);
T  = t1;
