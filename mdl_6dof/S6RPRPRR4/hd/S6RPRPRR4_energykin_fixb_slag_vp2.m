% Calculate kinetic energy for
% S6RPRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d5,d6,theta2]';
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
% Datum: 2019-03-09 03:46
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPRR4_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 03:44:20
% EndTime: 2019-03-09 03:44:21
% DurationCPUTime: 0.45s
% Computational Cost: add. (354->94), mult. (703->139), div. (0->0), fcn. (382->8), ass. (0->39)
t117 = m(3) / 0.2e1;
t116 = -pkin(3) - pkin(8);
t105 = sin(qJ(5));
t108 = cos(qJ(5));
t106 = sin(qJ(3));
t100 = t106 * qJD(1);
t109 = cos(qJ(3));
t102 = sin(pkin(10));
t94 = (pkin(1) * t102 + pkin(7)) * qJD(1);
t88 = qJD(2) * t109 - t106 * t94;
t113 = qJD(4) - t88;
t81 = pkin(4) * t100 + t116 * qJD(3) + t113;
t103 = cos(pkin(10));
t114 = -pkin(1) * t103 - pkin(2);
t112 = -qJ(4) * t106 + t114;
t84 = (t116 * t109 + t112) * qJD(1);
t75 = t105 * t81 + t108 * t84;
t89 = t106 * qJD(2) + t109 * t94;
t115 = qJD(1) * t109;
t97 = t100 + qJD(5);
t86 = -qJD(3) * qJ(4) - t89;
t74 = -t105 * t84 + t108 * t81;
t83 = pkin(4) * t115 - t86;
t107 = cos(qJ(6));
t104 = sin(qJ(6));
t96 = qJD(6) + t97;
t95 = t114 * qJD(1);
t93 = qJD(3) * t108 - t105 * t115;
t92 = -qJD(3) * t105 - t108 * t115;
t87 = (-pkin(3) * t109 + t112) * qJD(1);
t85 = -qJD(3) * pkin(3) + t113;
t80 = t104 * t92 + t107 * t93;
t79 = -t104 * t93 + t107 * t92;
t76 = -pkin(5) * t92 + t83;
t73 = pkin(9) * t92 + t75;
t72 = pkin(5) * t97 - pkin(9) * t93 + t74;
t71 = t104 * t72 + t107 * t73;
t70 = -t104 * t73 + t107 * t72;
t1 = qJD(2) ^ 2 * t117 + m(4) * (t88 ^ 2 + t89 ^ 2 + t95 ^ 2) / 0.2e1 + m(5) * (t85 ^ 2 + t86 ^ 2 + t87 ^ 2) / 0.2e1 + m(7) * (t70 ^ 2 + t71 ^ 2 + t76 ^ 2) / 0.2e1 + m(6) * (t74 ^ 2 + t75 ^ 2 + t83 ^ 2) / 0.2e1 + (t74 * mrSges(6,1) - t75 * mrSges(6,2) + Ifges(6,3) * t97 / 0.2e1) * t97 + (t70 * mrSges(7,1) - t71 * mrSges(7,2) + Ifges(7,3) * t96 / 0.2e1) * t96 + (t83 * mrSges(6,2) - t74 * mrSges(6,3) + Ifges(6,5) * t97 + Ifges(6,1) * t93 / 0.2e1) * t93 + (t76 * mrSges(7,2) - t70 * mrSges(7,3) + Ifges(7,5) * t96 + Ifges(7,1) * t80 / 0.2e1) * t80 + (-t83 * mrSges(6,1) + t75 * mrSges(6,3) + Ifges(6,4) * t93 + Ifges(6,6) * t97 + Ifges(6,2) * t92 / 0.2e1) * t92 + (-t76 * mrSges(7,1) + t71 * mrSges(7,3) + Ifges(7,4) * t80 + Ifges(7,6) * t96 + Ifges(7,2) * t79 / 0.2e1) * t79 + (t88 * mrSges(4,1) - t89 * mrSges(4,2) + t85 * mrSges(5,2) - t86 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3)) * qJD(3) + ((-t95 * mrSges(4,1) - t86 * mrSges(5,1) + t87 * mrSges(5,2) + t89 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t115) * t109 + (t85 * mrSges(5,1) + t95 * mrSges(4,2) - t88 * mrSges(4,3) - t87 * mrSges(5,3)) * t106 + ((-Ifges(5,5) + Ifges(4,6)) * t109 + (-Ifges(5,4) + Ifges(4,5)) * t106) * qJD(3) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t103 * mrSges(3,1) - t102 * mrSges(3,2) + (t102 ^ 2 + t103 ^ 2) * t117 * pkin(1)) * pkin(1) + ((Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t106 + (Ifges(4,4) + Ifges(5,6)) * t109) * t106) * qJD(1)) * qJD(1);
T  = t1;
