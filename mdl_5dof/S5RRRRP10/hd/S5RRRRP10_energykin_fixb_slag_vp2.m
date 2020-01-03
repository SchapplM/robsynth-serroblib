% Calculate kinetic energy for
% S5RRRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d4]';
% m_mdh [6x1]
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
% Datum: 2019-12-31 22:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRP10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:08:04
% EndTime: 2019-12-31 22:08:05
% DurationCPUTime: 0.34s
% Computational Cost: add. (490->86), mult. (1115->131), div. (0->0), fcn. (822->8), ass. (0->36)
t101 = sin(qJ(4));
t104 = cos(qJ(4));
t106 = cos(qJ(2));
t103 = sin(qJ(2));
t99 = sin(pkin(5));
t112 = t99 * qJD(1);
t109 = t103 * t112;
t111 = cos(pkin(5)) * qJD(1);
t110 = pkin(1) * t111;
t92 = -pkin(7) * t109 + t106 * t110;
t98 = qJD(2) + t111;
t86 = -pkin(2) * t98 - t92;
t102 = sin(qJ(3));
t105 = cos(qJ(3));
t90 = -t102 * t109 + t105 * t98;
t91 = t102 * t98 + t105 * t109;
t75 = -pkin(3) * t90 - pkin(9) * t91 + t86;
t108 = t106 * t112;
t93 = pkin(7) * t108 + t103 * t110;
t87 = pkin(8) * t98 + t93;
t88 = (-pkin(2) * t106 - pkin(8) * t103 - pkin(1)) * t112;
t80 = t102 * t88 + t105 * t87;
t94 = qJD(3) - t108;
t78 = pkin(9) * t94 + t80;
t71 = t101 * t75 + t104 * t78;
t70 = -t101 * t78 + t104 * t75;
t79 = -t102 * t87 + t105 * t88;
t77 = -pkin(3) * t94 - t79;
t107 = qJD(1) ^ 2;
t89 = qJD(4) - t90;
t82 = t101 * t94 + t104 * t91;
t81 = -t101 * t91 + t104 * t94;
t72 = -pkin(4) * t81 + qJD(5) + t77;
t69 = qJ(5) * t81 + t71;
t68 = pkin(4) * t89 - qJ(5) * t82 + t70;
t1 = m(4) * (t79 ^ 2 + t80 ^ 2 + t86 ^ 2) / 0.2e1 + t107 * Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 * t107 * t99 ^ 2 + t92 ^ 2 + t93 ^ 2) / 0.2e1 + m(6) * (t68 ^ 2 + t69 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + t71 ^ 2 + t77 ^ 2) / 0.2e1 + (t79 * mrSges(4,1) - t80 * mrSges(4,2) + Ifges(4,3) * t94 / 0.2e1) * t94 + (t86 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,5) * t94 + Ifges(4,1) * t91 / 0.2e1) * t91 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t106 / 0.2e1) * t106 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t106 + Ifges(3,1) * t103 / 0.2e1) * t103) * t112 + (-t92 * t103 + t93 * t106) * mrSges(3,3)) * t112 + (t92 * mrSges(3,1) - t93 * mrSges(3,2) + Ifges(3,3) * t98 / 0.2e1 + (Ifges(3,5) * t103 + Ifges(3,6) * t106) * t112) * t98 + (-t86 * mrSges(4,1) + t80 * mrSges(4,3) + Ifges(4,4) * t91 + Ifges(4,6) * t94 + Ifges(4,2) * t90 / 0.2e1) * t90 + (t70 * mrSges(5,1) + t68 * mrSges(6,1) - t71 * mrSges(5,2) - t69 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t89) * t89 + (t77 * mrSges(5,2) + t72 * mrSges(6,2) - t70 * mrSges(5,3) - t68 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t82 + (Ifges(5,5) + Ifges(6,5)) * t89) * t82 + (-t77 * mrSges(5,1) - t72 * mrSges(6,1) + t71 * mrSges(5,3) + t69 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t81 + (Ifges(5,6) + Ifges(6,6)) * t89 + (Ifges(5,4) + Ifges(6,4)) * t82) * t81;
T = t1;
