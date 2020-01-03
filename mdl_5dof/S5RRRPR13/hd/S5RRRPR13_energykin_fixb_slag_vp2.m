% Calculate kinetic energy for
% S5RRRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d1,d2,d3,d5]';
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
% Datum: 2019-12-31 21:48
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR13_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 21:43:13
% EndTime: 2019-12-31 21:43:14
% DurationCPUTime: 0.36s
% Computational Cost: add. (420->87), mult. (959->131), div. (0->0), fcn. (684->8), ass. (0->39)
t114 = pkin(3) + pkin(9);
t113 = cos(qJ(3));
t101 = sin(qJ(3));
t102 = sin(qJ(2));
t104 = cos(qJ(2));
t98 = sin(pkin(5));
t111 = t98 * qJD(1);
t108 = t104 * t111;
t112 = cos(pkin(5)) * qJD(1);
t110 = pkin(1) * t112;
t91 = pkin(7) * t108 + t102 * t110;
t97 = qJD(2) + t112;
t85 = pkin(8) * t97 + t91;
t86 = (-pkin(2) * t104 - pkin(8) * t102 - pkin(1)) * t111;
t78 = t101 * t86 + t113 * t85;
t109 = t102 * t111;
t93 = qJD(3) - t108;
t76 = -qJ(4) * t93 - t78;
t77 = -t101 * t85 + t113 * t86;
t90 = -pkin(7) * t109 + t104 * t110;
t107 = qJD(4) - t77;
t84 = -pkin(2) * t97 - t90;
t89 = t101 * t97 + t113 * t109;
t106 = -qJ(4) * t89 + t84;
t105 = qJD(1) ^ 2;
t103 = cos(qJ(5));
t100 = sin(qJ(5));
t88 = t101 * t109 - t113 * t97;
t87 = qJD(5) + t89;
t80 = t100 * t88 + t103 * t93;
t79 = -t100 * t93 + t103 * t88;
t75 = -t93 * pkin(3) + t107;
t74 = pkin(3) * t88 + t106;
t73 = -pkin(4) * t88 - t76;
t72 = t114 * t88 + t106;
t71 = t89 * pkin(4) - t114 * t93 + t107;
t70 = t100 * t71 + t103 * t72;
t69 = -t100 * t72 + t103 * t71;
t1 = m(3) * (pkin(1) ^ 2 * t105 * t98 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + t105 * Ifges(2,3) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + m(4) * (t77 ^ 2 + t78 ^ 2 + t84 ^ 2) / 0.2e1 + m(6) * (t69 ^ 2 + t70 ^ 2 + t73 ^ 2) / 0.2e1 + (t69 * mrSges(6,1) - t70 * mrSges(6,2) + Ifges(6,3) * t87 / 0.2e1) * t87 + (t73 * mrSges(6,2) - t69 * mrSges(6,3) + Ifges(6,5) * t87 + Ifges(6,1) * t80 / 0.2e1) * t80 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t104 / 0.2e1) * t104 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t104 + Ifges(3,1) * t102 / 0.2e1) * t102) * t111 + (-t90 * t102 + t91 * t104) * mrSges(3,3)) * t111 + (t90 * mrSges(3,1) - t91 * mrSges(3,2) + Ifges(3,3) * t97 / 0.2e1 + (Ifges(3,5) * t102 + Ifges(3,6) * t104) * t111) * t97 + (-t73 * mrSges(6,1) + t70 * mrSges(6,3) + Ifges(6,4) * t80 + Ifges(6,6) * t87 + Ifges(6,2) * t79 / 0.2e1) * t79 + (t77 * mrSges(4,1) - t78 * mrSges(4,2) + t75 * mrSges(5,2) - t76 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,3) / 0.2e1) * t93) * t93 + (t75 * mrSges(5,1) + t84 * mrSges(4,2) - t77 * mrSges(4,3) - t74 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1) * t89 + (-Ifges(5,4) + Ifges(4,5)) * t93) * t89 + (t84 * mrSges(4,1) + t76 * mrSges(5,1) - t74 * mrSges(5,2) - t78 * mrSges(4,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * t88 + (Ifges(5,5) - Ifges(4,6)) * t93 + (-Ifges(4,4) - Ifges(5,6)) * t89) * t88;
T = t1;
