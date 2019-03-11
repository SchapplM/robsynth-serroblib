% Calculate kinetic energy for
% S6RPRRRP1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 05:57
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRRP1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP1_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:55:22
% EndTime: 2019-03-09 05:55:22
% DurationCPUTime: 0.52s
% Computational Cost: add. (444->94), mult. (899->140), div. (0->0), fcn. (570->8), ass. (0->36)
t104 = sin(qJ(4));
t105 = sin(qJ(3));
t106 = cos(qJ(4));
t107 = cos(qJ(3));
t91 = (t104 * t105 - t106 * t107) * qJD(1);
t114 = m(3) / 0.2e1;
t113 = cos(qJ(5));
t103 = sin(qJ(5));
t100 = qJD(3) + qJD(4);
t112 = pkin(8) * qJD(1);
t101 = sin(pkin(10));
t95 = (pkin(1) * t101 + pkin(7)) * qJD(1);
t99 = t107 * qJD(2);
t84 = qJD(3) * pkin(3) + t99 + (-t95 - t112) * t105;
t89 = t105 * qJD(2) + t107 * t95;
t87 = t107 * t112 + t89;
t78 = t104 * t84 + t106 * t87;
t76 = pkin(9) * t100 + t78;
t92 = (t104 * t107 + t105 * t106) * qJD(1);
t102 = cos(pkin(10));
t111 = -pkin(1) * t102 - pkin(2);
t93 = (-pkin(3) * t107 + t111) * qJD(1);
t80 = pkin(4) * t91 - pkin(9) * t92 + t93;
t72 = t103 * t80 + t113 * t76;
t77 = -t104 * t87 + t106 * t84;
t75 = -pkin(4) * t100 - t77;
t71 = -t103 * t76 + t113 * t80;
t96 = t111 * qJD(1);
t90 = qJD(5) + t91;
t88 = -t105 * t95 + t99;
t86 = t103 * t100 + t113 * t92;
t85 = -t100 * t113 + t103 * t92;
t73 = pkin(5) * t85 - qJ(6) * t86 + t75;
t70 = qJ(6) * t90 + t72;
t69 = -t90 * pkin(5) + qJD(6) - t71;
t1 = qJD(2) ^ 2 * t114 + m(5) * (t77 ^ 2 + t78 ^ 2 + t93 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t96 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + m(6) * (t71 ^ 2 + t72 ^ 2 + t75 ^ 2) / 0.2e1 + (t93 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 - (-t93 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t92 - Ifges(5,2) * t91 / 0.2e1) * t91 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + Ifges(5,5) * t92 - Ifges(5,6) * t91 + Ifges(5,3) * t100 / 0.2e1) * t100 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t71 * mrSges(6,1) - t69 * mrSges(7,1) - t72 * mrSges(6,2) + t70 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t90) * t90 + (t75 * mrSges(6,2) + t69 * mrSges(7,2) - t71 * mrSges(6,3) - t73 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t86 + (Ifges(7,4) + Ifges(6,5)) * t90) * t86 + (t75 * mrSges(6,1) + t73 * mrSges(7,1) - t70 * mrSges(7,2) - t72 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t85 + (-Ifges(6,6) + Ifges(7,6)) * t90 + (-Ifges(6,4) + Ifges(7,5)) * t86) * t85 + (t96 * (-mrSges(4,1) * t107 + mrSges(4,2) * t105) + (-t88 * t105 + t89 * t107) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t105 + Ifges(4,6) * t107) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t102 * mrSges(3,1) - t101 * mrSges(3,2) + (t101 ^ 2 + t102 ^ 2) * t114 * pkin(1)) * pkin(1) + Ifges(4,2) * t107 ^ 2 / 0.2e1 + (Ifges(4,4) * t107 + Ifges(4,1) * t105 / 0.2e1) * t105) * qJD(1)) * qJD(1);
T  = t1;
