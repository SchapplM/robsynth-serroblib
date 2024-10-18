% Calculate kinetic energy for
% S5PRRRR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha3,d2,d3,d4,d5,theta1]';
% m [6x1]
%   mass of all robot links (including the base)
% mrSges [6x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [6x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2024-09-27 21:46
% Revision: 4128c0dfad36ceb08eaa7a9fbc8f797c888f59be (2024-07-29)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR11_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR11_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2024-09-27 21:45:08
% EndTime: 2024-09-27 21:45:08
% DurationCPUTime: 0.06s
% Computational Cost: add. (302->72), mult. (838->120), div. (0->0), fcn. (626->8), ass. (0->36)
t105 = cos(pkin(5));
t102 = t105 * qJD(2);
t104 = sin(pkin(5));
t119 = pkin(2) * t102 + qJD(1) * t104;
t107 = sin(qJ(4));
t110 = cos(qJ(4));
t101 = t102 + qJD(3);
t108 = sin(qJ(3));
t116 = t104 * qJD(2);
t114 = t108 * t116;
t111 = cos(qJ(3));
t118 = t119 * t111;
t85 = t101 * pkin(3) + (-pkin(7) - pkin(8)) * t114 + t118;
t113 = t111 * t116;
t89 = pkin(7) * t113 + t119 * t108;
t87 = pkin(8) * t113 + t89;
t79 = t107 * t85 + t110 * t87;
t95 = qJD(4) + t101;
t78 = -t107 * t87 + t110 * t85;
t103 = t105 * qJD(1);
t92 = t103 + (-pkin(3) * t111 - pkin(2)) * t116;
t109 = cos(qJ(5));
t106 = sin(qJ(5));
t94 = qJD(5) + t95;
t93 = -pkin(2) * t116 + t103;
t91 = (t107 * t111 + t108 * t110) * t116;
t90 = (-t107 * t108 + t110 * t111) * t116;
t88 = -pkin(7) * t114 + t118;
t84 = -t90 * pkin(4) + t92;
t81 = t106 * t90 + t109 * t91;
t80 = -t106 * t91 + t109 * t90;
t77 = t90 * pkin(9) + t79;
t76 = t95 * pkin(4) - t91 * pkin(9) + t78;
t75 = t106 * t76 + t109 * t77;
t74 = -t106 * t77 + t109 * t76;
t1 = m(5) * (t78 ^ 2 + t79 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t84 ^ 2) / 0.2e1 + m(4) * (t88 ^ 2 + t89 ^ 2 + t93 ^ 2) / 0.2e1 + (m(3) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t78 * mrSges(5,1) - t79 * mrSges(5,2) + Ifges(5,3) * t95 / 0.2e1) * t95 + (t74 * mrSges(6,1) - t75 * mrSges(6,2) + Ifges(6,3) * t94 / 0.2e1) * t94 + (t92 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,5) * t95 + Ifges(5,1) * t91 / 0.2e1) * t91 + (t84 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,5) * t94 + Ifges(6,1) * t81 / 0.2e1) * t81 + (-t92 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t91 + Ifges(5,6) * t95 + Ifges(5,2) * t90 / 0.2e1) * t90 + (-t84 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t81 + Ifges(6,6) * t94 + Ifges(6,2) * t80 / 0.2e1) * t80 + (Ifges(3,3) * qJD(2) / 0.2e1 + (t93 * (-mrSges(4,1) * t111 + mrSges(4,2) * t108) + (Ifges(4,2) * t111 ^ 2 / 0.2e1 + (Ifges(4,4) * t111 + Ifges(4,1) * t108 / 0.2e1) * t108) * t116 + (-t88 * t108 + t89 * t111) * mrSges(4,3)) * t104) * qJD(2) + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + Ifges(4,3) * t101 / 0.2e1 + (Ifges(4,5) * t108 + Ifges(4,6) * t111) * t116) * t101;
T = t1;
