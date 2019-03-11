% Calculate kinetic energy for
% S6RPRPRP5
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
% Datum: 2019-03-09 03:17
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:14:22
% EndTime: 2019-03-09 03:14:22
% DurationCPUTime: 0.55s
% Computational Cost: add. (659->100), mult. (1599->145), div. (0->0), fcn. (1182->8), ass. (0->40)
t110 = cos(pkin(9));
t122 = t110 ^ 2;
t121 = m(3) / 0.2e1;
t120 = cos(qJ(5));
t119 = pkin(7) + qJ(2);
t111 = sin(qJ(5));
t107 = sin(pkin(10));
t109 = cos(pkin(10));
t102 = qJD(2) + (-pkin(2) * t110 - pkin(1)) * qJD(1);
t112 = sin(qJ(3));
t113 = cos(qJ(3));
t116 = t110 * qJD(1);
t108 = sin(pkin(9));
t117 = t108 * qJD(1);
t98 = t112 * t117 - t113 * t116;
t99 = (t108 * t113 + t110 * t112) * qJD(1);
t85 = pkin(3) * t98 - qJ(4) * t99 + t102;
t100 = t119 * t117;
t101 = t119 * t116;
t90 = -t112 * t100 + t113 * t101;
t88 = qJD(3) * qJ(4) + t90;
t78 = -t107 * t88 + t109 * t85;
t93 = qJD(3) * t107 + t109 * t99;
t75 = pkin(4) * t98 - pkin(8) * t93 + t78;
t79 = t107 * t85 + t109 * t88;
t92 = qJD(3) * t109 - t107 * t99;
t77 = pkin(8) * t92 + t79;
t72 = t111 * t75 + t120 * t77;
t89 = -t100 * t113 - t112 * t101;
t71 = -t111 * t77 + t120 * t75;
t87 = -qJD(3) * pkin(3) + qJD(4) - t89;
t80 = -pkin(4) * t92 + t87;
t104 = -qJD(1) * pkin(1) + qJD(2);
t94 = qJD(5) + t98;
t82 = t111 * t92 + t120 * t93;
t81 = t111 * t93 - t120 * t92;
t73 = pkin(5) * t81 - qJ(6) * t82 + t80;
t70 = qJ(6) * t94 + t72;
t69 = -t94 * pkin(5) + qJD(6) - t71;
t1 = t104 ^ 2 * t121 + m(4) * (t102 ^ 2 + t89 ^ 2 + t90 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t80 ^ 2) / 0.2e1 + (t102 * mrSges(4,2) - t89 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t99 + (t87 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (-t87 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (t89 * mrSges(4,1) - t90 * mrSges(4,2) + Ifges(4,5) * t99 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t94) * t94 + (t80 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t82 + (Ifges(7,4) + Ifges(6,5)) * t94) * t82 + (t102 * mrSges(4,1) + t78 * mrSges(5,1) - t79 * mrSges(5,2) - t90 * mrSges(4,3) - Ifges(4,4) * t99 + Ifges(5,5) * t93 - Ifges(4,6) * qJD(3) + Ifges(5,6) * t92 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t98) * t98 + (t80 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t81 + (-Ifges(6,6) + Ifges(7,6)) * t94 + (-Ifges(6,4) + Ifges(7,5)) * t82) * t81 + (t104 * (-mrSges(3,1) * t110 + mrSges(3,2) * t108) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t121 + mrSges(3,3)) * (t108 ^ 2 + t122) * qJ(2) + Ifges(3,2) * t122 / 0.2e1 + (Ifges(3,4) * t110 + Ifges(3,1) * t108 / 0.2e1) * t108) * qJD(1)) * qJD(1);
T  = t1;
