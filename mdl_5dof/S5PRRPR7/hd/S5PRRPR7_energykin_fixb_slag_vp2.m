% Calculate kinetic energy for
% S5PRRPR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,d2,d3,d5,theta1,theta4]';
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
% Datum: 2019-12-05 16:38
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRPR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(10,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRPR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 16:34:54
% EndTime: 2019-12-05 16:34:55
% DurationCPUTime: 0.29s
% Computational Cost: add. (238->74), mult. (544->121), div. (0->0), fcn. (355->10), ass. (0->35)
t100 = cos(pkin(10));
t103 = sin(qJ(3));
t106 = cos(qJ(3));
t101 = cos(pkin(5));
t112 = qJD(1) * t101;
t104 = sin(qJ(2));
t99 = sin(pkin(5));
t113 = qJD(1) * t99;
t94 = qJD(2) * pkin(7) + t104 * t113;
t87 = t103 * t112 + t106 * t94;
t83 = qJD(3) * qJ(4) + t87;
t107 = cos(qJ(2));
t109 = t107 * t113;
t88 = -t109 + (-pkin(3) * t106 - qJ(4) * t103 - pkin(2)) * qJD(2);
t98 = sin(pkin(10));
t79 = t100 * t83 + t98 * t88;
t111 = qJD(2) * t103;
t110 = qJD(2) * t106;
t86 = -t103 * t94 + t106 * t112;
t78 = t100 * t88 - t98 * t83;
t92 = t100 * qJD(3) - t98 * t111;
t81 = -qJD(3) * pkin(3) + qJD(4) - t86;
t105 = cos(qJ(5));
t102 = sin(qJ(5));
t95 = -qJD(2) * pkin(2) - t109;
t93 = t98 * qJD(3) + t100 * t111;
t89 = qJD(5) - t92;
t85 = -t102 * t110 + t105 * t93;
t84 = -t102 * t93 - t105 * t110;
t77 = -t92 * pkin(4) - t93 * pkin(8) + t81;
t76 = -pkin(8) * t110 + t79;
t75 = pkin(4) * t110 - t78;
t74 = t102 * t77 + t105 * t76;
t73 = -t102 * t76 + t105 * t77;
t1 = m(4) * (t86 ^ 2 + t87 ^ 2 + t95 ^ 2) / 0.2e1 + m(6) * (t73 ^ 2 + t74 ^ 2 + t75 ^ 2) / 0.2e1 + m(5) * (t78 ^ 2 + t79 ^ 2 + t81 ^ 2) / 0.2e1 + (t81 * mrSges(5,2) - t78 * mrSges(5,3) + Ifges(5,1) * t93 / 0.2e1) * t93 + (t73 * mrSges(6,1) - t74 * mrSges(6,2) + Ifges(6,3) * t89 / 0.2e1) * t89 + (t86 * mrSges(4,1) - t87 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (-t81 * mrSges(5,1) + t79 * mrSges(5,3) + Ifges(5,4) * t93 + Ifges(5,2) * t92 / 0.2e1) * t92 + (t75 * mrSges(6,2) - t73 * mrSges(6,3) + Ifges(6,5) * t89 + Ifges(6,1) * t85 / 0.2e1) * t85 + (-t75 * mrSges(6,1) + t74 * mrSges(6,3) + Ifges(6,4) * t85 + Ifges(6,6) * t89 + Ifges(6,2) * t84 / 0.2e1) * t84 + (m(3) * (t101 ^ 2 + (t104 ^ 2 + t107 ^ 2) * t99 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * qJD(2) / 0.2e1 + (mrSges(3,1) * t107 - mrSges(3,2) * t104) * t113 + (t95 * mrSges(4,2) - t86 * mrSges(4,3) + Ifges(4,5) * qJD(3) + Ifges(4,1) * t111 / 0.2e1) * t103 + (-t95 * mrSges(4,1) - t78 * mrSges(5,1) + t79 * mrSges(5,2) + t87 * mrSges(4,3) - Ifges(5,5) * t93 + Ifges(4,6) * qJD(3) - Ifges(5,6) * t92 + (Ifges(4,4) * t103 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t106) * qJD(2)) * t106) * qJD(2);
T = t1;
