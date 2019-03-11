% Calculate kinetic energy for
% S6RPPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d4,d6,theta2,theta3]';
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
% Datum: 2019-03-09 01:42
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPPRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPPRPR2_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 01:41:28
% EndTime: 2019-03-09 01:41:29
% DurationCPUTime: 0.42s
% Computational Cost: add. (333->90), mult. (758->132), div. (0->0), fcn. (472->8), ass. (0->38)
t113 = m(3) / 0.2e1;
t112 = pkin(4) + pkin(8);
t111 = cos(qJ(4));
t103 = sin(qJ(4));
t99 = sin(pkin(9));
t94 = (pkin(1) * t99 + qJ(3)) * qJD(1);
t100 = cos(pkin(10));
t97 = t100 * qJD(2);
t98 = sin(pkin(10));
t82 = t97 + (-pkin(7) * qJD(1) - t94) * t98;
t110 = qJD(1) * t100;
t87 = t98 * qJD(2) + t100 * t94;
t83 = pkin(7) * t110 + t87;
t77 = t103 * t82 + t111 * t83;
t101 = cos(pkin(9));
t109 = -pkin(1) * t101 - pkin(2);
t76 = -t103 * t83 + t111 * t82;
t75 = -qJD(4) * qJ(5) - t77;
t108 = qJD(5) - t76;
t91 = (t100 * t103 + t111 * t98) * qJD(1);
t89 = qJD(3) + (-pkin(3) * t100 + t109) * qJD(1);
t107 = -qJ(5) * t91 + t89;
t104 = cos(qJ(6));
t102 = sin(qJ(6));
t93 = qJD(1) * t109 + qJD(3);
t90 = qJD(1) * t103 * t98 - t111 * t110;
t88 = qJD(6) + t91;
t86 = -t94 * t98 + t97;
t85 = qJD(4) * t104 + t102 * t90;
t84 = -qJD(4) * t102 + t104 * t90;
t78 = pkin(4) * t90 + t107;
t74 = -qJD(4) * pkin(4) + t108;
t73 = t112 * t90 + t107;
t72 = -pkin(5) * t90 - t75;
t71 = t91 * pkin(5) - t112 * qJD(4) + t108;
t70 = t102 * t71 + t104 * t73;
t69 = -t102 * t73 + t104 * t71;
t1 = qJD(2) ^ 2 * t113 + m(4) * (t86 ^ 2 + t87 ^ 2 + t93 ^ 2) / 0.2e1 + m(5) * (t76 ^ 2 + t77 ^ 2 + t89 ^ 2) / 0.2e1 + m(7) * (t69 ^ 2 + t70 ^ 2 + t72 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t78 ^ 2) / 0.2e1 + (t69 * mrSges(7,1) - t70 * mrSges(7,2) + Ifges(7,3) * t88 / 0.2e1) * t88 + (t72 * mrSges(7,2) - t69 * mrSges(7,3) + Ifges(7,5) * t88 + Ifges(7,1) * t85 / 0.2e1) * t85 + (-t72 * mrSges(7,1) + t70 * mrSges(7,3) + Ifges(7,4) * t85 + Ifges(7,6) * t88 + Ifges(7,2) * t84 / 0.2e1) * t84 + (t74 * mrSges(6,1) + t89 * mrSges(5,2) - t76 * mrSges(5,3) - t78 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t91) * t91 + (t89 * mrSges(5,1) + t75 * mrSges(6,1) - t78 * mrSges(6,2) - t77 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,3) / 0.2e1) * t90 + (-Ifges(5,4) - Ifges(6,6)) * t91) * t90 + (t76 * mrSges(5,1) - t77 * mrSges(5,2) + t74 * mrSges(6,2) - t75 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(6,1) / 0.2e1) * qJD(4) + (-Ifges(6,4) + Ifges(5,5)) * t91 + (Ifges(6,5) - Ifges(5,6)) * t90) * qJD(4) + (t93 * (-mrSges(4,1) * t100 + mrSges(4,2) * t98) + (t87 * t100 - t86 * t98) * mrSges(4,3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t101 * mrSges(3,1) - t99 * mrSges(3,2) + (t101 ^ 2 + t99 ^ 2) * t113 * pkin(1)) * pkin(1) + Ifges(4,1) * t98 ^ 2 / 0.2e1 + (Ifges(4,4) * t98 + Ifges(4,2) * t100 / 0.2e1) * t100) * qJD(1)) * qJD(1);
T  = t1;
