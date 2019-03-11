% Calculate kinetic energy for
% S6RPRPRP3
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
% Datum: 2019-03-09 03:10
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRP3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRP3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:07:44
% EndTime: 2019-03-09 03:07:44
% DurationCPUTime: 0.47s
% Computational Cost: add. (430->93), mult. (919->136), div. (0->0), fcn. (574->8), ass. (0->35)
t112 = m(3) / 0.2e1;
t111 = cos(qJ(5));
t103 = sin(qJ(5));
t105 = cos(qJ(3));
t109 = t105 * qJD(1);
t101 = cos(pkin(10));
t104 = sin(qJ(3));
t100 = sin(pkin(9));
t95 = (pkin(1) * t100 + pkin(7)) * qJD(1);
t89 = t104 * qJD(2) + t105 * t95;
t86 = qJD(3) * qJ(4) + t89;
t102 = cos(pkin(9));
t108 = -pkin(1) * t102 - pkin(2);
t87 = (-pkin(3) * t105 - qJ(4) * t104 + t108) * qJD(1);
t99 = sin(pkin(10));
t77 = t101 * t87 - t86 * t99;
t110 = t104 * qJD(1);
t92 = qJD(3) * t99 + t101 * t110;
t74 = -pkin(4) * t109 - pkin(8) * t92 + t77;
t78 = t101 * t86 + t99 * t87;
t91 = qJD(3) * t101 - t110 * t99;
t76 = pkin(8) * t91 + t78;
t71 = t103 * t74 + t111 * t76;
t88 = qJD(2) * t105 - t104 * t95;
t70 = -t103 * t76 + t111 * t74;
t85 = -qJD(3) * pkin(3) + qJD(4) - t88;
t79 = -pkin(4) * t91 + t85;
t97 = qJD(5) - t109;
t96 = t108 * qJD(1);
t81 = t103 * t91 + t111 * t92;
t80 = t103 * t92 - t111 * t91;
t72 = pkin(5) * t80 - qJ(6) * t81 + t79;
t69 = qJ(6) * t97 + t71;
t68 = -t97 * pkin(5) + qJD(6) - t70;
t1 = qJD(2) ^ 2 * t112 + m(4) * (t88 ^ 2 + t89 ^ 2 + t96 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t71 ^ 2 + t79 ^ 2) / 0.2e1 + m(7) * (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) / 0.2e1 + (t85 * mrSges(5,2) - t77 * mrSges(5,3) + Ifges(5,1) * t92 / 0.2e1) * t92 + (-t85 * mrSges(5,1) + t78 * mrSges(5,3) + Ifges(5,4) * t92 + Ifges(5,2) * t91 / 0.2e1) * t91 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t70 * mrSges(6,1) - t68 * mrSges(7,1) - t71 * mrSges(6,2) + t69 * mrSges(7,3) + (Ifges(7,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t97) * t97 + (t79 * mrSges(6,2) + t68 * mrSges(7,2) - t70 * mrSges(6,3) - t72 * mrSges(7,3) + (Ifges(7,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t81 + (Ifges(7,4) + Ifges(6,5)) * t97) * t81 + (t79 * mrSges(6,1) + t72 * mrSges(7,1) - t69 * mrSges(7,2) - t71 * mrSges(6,3) + (Ifges(7,3) / 0.2e1 + Ifges(6,2) / 0.2e1) * t80 + (-Ifges(6,6) + Ifges(7,6)) * t97 + (-Ifges(6,4) + Ifges(7,5)) * t81) * t80 + ((t96 * mrSges(4,2) - t88 * mrSges(4,3) + Ifges(4,1) * t110 / 0.2e1) * t104 + (-t96 * mrSges(4,1) - t77 * mrSges(5,1) + t78 * mrSges(5,2) + t89 * mrSges(4,3) - Ifges(5,5) * t92 - Ifges(5,6) * t91) * t105 + qJD(3) * (Ifges(4,5) * t104 + Ifges(4,6) * t105) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t102 * mrSges(3,1) - t100 * mrSges(3,2) + (t100 ^ 2 + t102 ^ 2) * t112 * pkin(1)) * pkin(1) + (Ifges(4,4) * t104 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t105) * t105) * qJD(1)) * qJD(1);
T  = t1;
