% Calculate kinetic energy for
% S6RPRRPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,d6,theta2]';
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
% Datum: 2019-03-09 05:14
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRRPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 05:12:07
% EndTime: 2019-03-09 05:12:07
% DurationCPUTime: 0.58s
% Computational Cost: add. (607->101), mult. (1491->147), div. (0->0), fcn. (1106->8), ass. (0->42)
t107 = cos(pkin(10));
t123 = t107 ^ 2;
t122 = qJD(1) * (pkin(7) + qJ(2));
t121 = m(3) / 0.2e1;
t120 = pkin(4) + pkin(9);
t119 = cos(qJ(4));
t109 = sin(qJ(4));
t100 = t107 * t122;
t110 = sin(qJ(3));
t112 = cos(qJ(3));
t106 = sin(pkin(10));
t99 = t106 * t122;
t90 = -t100 * t110 - t112 * t99;
t98 = (t106 * t112 + t107 * t110) * qJD(1);
t83 = qJD(3) * pkin(3) - pkin(8) * t98 + t90;
t91 = t112 * t100 - t110 * t99;
t97 = (-t106 * t110 + t107 * t112) * qJD(1);
t84 = pkin(8) * t97 + t91;
t78 = t109 * t83 + t119 * t84;
t105 = qJD(3) + qJD(4);
t76 = -qJ(5) * t105 - t78;
t77 = -t109 * t84 + t119 * t83;
t116 = qJD(5) - t77;
t89 = t109 * t97 + t119 * t98;
t101 = qJD(2) + (-pkin(2) * t107 - pkin(1)) * qJD(1);
t92 = -pkin(3) * t97 + t101;
t115 = -qJ(5) * t89 + t92;
t111 = cos(qJ(6));
t108 = sin(qJ(6));
t102 = -qJD(1) * pkin(1) + qJD(2);
t88 = t109 * t98 - t119 * t97;
t87 = qJD(6) + t89;
t86 = t105 * t111 + t108 * t88;
t85 = -t105 * t108 + t111 * t88;
t79 = pkin(4) * t88 + t115;
t75 = -t105 * pkin(4) + t116;
t74 = t120 * t88 + t115;
t73 = -pkin(5) * t88 - t76;
t72 = t89 * pkin(5) - t120 * t105 + t116;
t71 = t108 * t72 + t111 * t74;
t70 = -t108 * t74 + t111 * t72;
t1 = t102 ^ 2 * t121 + m(4) * (t101 ^ 2 + t90 ^ 2 + t91 ^ 2) / 0.2e1 + m(5) * (t77 ^ 2 + t78 ^ 2 + t92 ^ 2) / 0.2e1 + m(6) * (t75 ^ 2 + t76 ^ 2 + t79 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t73 ^ 2) / 0.2e1 + (t101 * mrSges(4,2) - t90 * mrSges(4,3) + Ifges(4,1) * t98 / 0.2e1) * t98 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t87 / 0.2e1) * t87 + (-t101 * mrSges(4,1) + t91 * mrSges(4,3) + Ifges(4,4) * t98 + Ifges(4,2) * t97 / 0.2e1) * t97 + (t73 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t87 + Ifges(7,1) * t86 / 0.2e1) * t86 + (-t73 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t86 + Ifges(7,6) * t87 + Ifges(7,2) * t85 / 0.2e1) * t85 + (t90 * mrSges(4,1) - t91 * mrSges(4,2) + Ifges(4,5) * t98 + Ifges(4,6) * t97 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t75 * mrSges(6,1) + t92 * mrSges(5,2) - t77 * mrSges(5,3) - t79 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1) * t89) * t89 + (t92 * mrSges(5,1) + t76 * mrSges(6,1) - t79 * mrSges(6,2) - t78 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1) * t88 + (-Ifges(5,4) - Ifges(6,6)) * t89) * t88 + (t77 * mrSges(5,1) - t78 * mrSges(5,2) + t75 * mrSges(6,2) - t76 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1) * t105 + (-Ifges(6,4) + Ifges(5,5)) * t89 + (Ifges(6,5) - Ifges(5,6)) * t88) * t105 + (t102 * (-mrSges(3,1) * t107 + mrSges(3,2) * t106) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t121 + mrSges(3,3)) * (t106 ^ 2 + t123) * qJ(2) + Ifges(3,2) * t123 / 0.2e1 + (Ifges(3,4) * t107 + Ifges(3,1) * t106 / 0.2e1) * t106) * qJD(1)) * qJD(1);
T  = t1;
