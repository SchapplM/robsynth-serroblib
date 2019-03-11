% Calculate kinetic energy for
% S6RRPPPR5
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d2,d6,theta3]';
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
% Datum: 2019-03-09 08:24
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RRPPPR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RRPPPR5_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 08:21:16
% EndTime: 2019-03-09 08:21:16
% DurationCPUTime: 0.49s
% Computational Cost: add. (384->108), mult. (839->138), div. (0->0), fcn. (500->6), ass. (0->36)
t97 = sin(qJ(2));
t106 = t97 * qJD(1);
t107 = cos(pkin(9));
t95 = sin(pkin(9));
t84 = -t107 * qJD(2) + t95 * t106;
t112 = (pkin(3) + qJ(5)) * t84;
t111 = pkin(4) + pkin(8);
t110 = pkin(7) * mrSges(3,3);
t108 = -pkin(5) - qJ(4);
t99 = cos(qJ(2));
t83 = (-pkin(2) * t99 - qJ(3) * t97 - pkin(1)) * qJD(1);
t105 = t99 * qJD(1);
t89 = pkin(7) * t105 + qJD(2) * qJ(3);
t81 = t107 * t89 + t95 * t83;
t104 = qJD(2) * pkin(2) - pkin(7) * t106 - qJD(3);
t80 = t107 * t83 - t95 * t89;
t77 = qJ(4) * t105 - t81;
t76 = pkin(3) * t105 + qJD(4) - t80;
t85 = t95 * qJD(2) + t107 * t106;
t103 = qJ(4) * t85 + t104;
t102 = qJ(5) * t105 + t76;
t98 = cos(qJ(6));
t96 = sin(qJ(6));
t90 = qJD(6) - t105;
t79 = t84 * t98 + t85 * t96;
t78 = -t84 * t96 + t85 * t98;
t75 = pkin(3) * t84 - t103;
t74 = -pkin(4) * t84 + qJD(5) - t77;
t73 = t103 - t112;
t72 = t85 * pkin(4) + t102;
t71 = t108 * t85 - t104 + t112;
t70 = t108 * t105 - t111 * t84 + qJD(5) + t81;
t69 = t111 * t85 + t102;
t68 = t69 * t98 + t70 * t96;
t67 = -t69 * t96 + t70 * t98;
t1 = m(4) * (t104 ^ 2 + t80 ^ 2 + t81 ^ 2) / 0.2e1 + Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(7) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t72 ^ 2 + t73 ^ 2 + t74 ^ 2) / 0.2e1 + m(5) * (t75 ^ 2 + t76 ^ 2 + t77 ^ 2) / 0.2e1 + (t67 * mrSges(7,1) - t68 * mrSges(7,2) + Ifges(7,3) * t90 / 0.2e1) * t90 + (t71 * mrSges(7,2) - t67 * mrSges(7,3) + Ifges(7,5) * t90 + Ifges(7,1) * t79 / 0.2e1) * t79 + (-t71 * mrSges(7,1) + t68 * mrSges(7,3) + Ifges(7,4) * t79 + Ifges(7,6) * t90 + Ifges(7,2) * t78 / 0.2e1) * t78 + (t76 * mrSges(5,1) + t73 * mrSges(6,1) - t104 * mrSges(4,2) - t72 * mrSges(6,2) - t80 * mrSges(4,3) - t75 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,3) / 0.2e1) * t85) * t85 + (-t104 * mrSges(4,1) + t77 * mrSges(5,1) - t75 * mrSges(5,2) + t74 * mrSges(6,2) - t81 * mrSges(4,3) - t73 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t84 + (-Ifges(4,4) + Ifges(6,5) - Ifges(5,6)) * t85 + (Ifges(6,4) - Ifges(5,5) + Ifges(4,6)) * t105) * t84 + ((-pkin(7) * mrSges(3,1) + Ifges(3,5)) * qJD(2) * t97 + (-t80 * mrSges(4,1) - t74 * mrSges(6,1) + t81 * mrSges(4,2) - t76 * mrSges(5,2) + t77 * mrSges(5,3) + t72 * mrSges(6,3) + (-pkin(7) * mrSges(3,2) + Ifges(3,6)) * qJD(2) + (Ifges(5,4) + Ifges(6,6) - Ifges(4,5)) * t85) * t99 + (Ifges(2,3) / 0.2e1 + m(3) * (pkin(1) ^ 2 + (t97 ^ 2 + t99 ^ 2) * pkin(7) ^ 2) / 0.2e1 + (-pkin(1) * mrSges(3,2) + (t110 + Ifges(3,1) / 0.2e1) * t97) * t97 + (pkin(1) * mrSges(3,1) + Ifges(3,4) * t97 + (t110 + Ifges(3,2) / 0.2e1 + Ifges(6,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * t99) * t99) * qJD(1)) * qJD(1);
T  = t1;
