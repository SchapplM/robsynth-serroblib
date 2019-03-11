% Calculate kinetic energy for
% S6RPRPPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6,theta2]';
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
% Datum: 2019-03-09 02:45
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(9,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR3_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:44:02
% EndTime: 2019-03-09 02:44:02
% DurationCPUTime: 0.37s
% Computational Cost: add. (226->91), mult. (455->121), div. (0->0), fcn. (194->6), ass. (0->32)
t107 = m(3) / 0.2e1;
t106 = -pkin(3) - pkin(4);
t92 = sin(pkin(9));
t83 = (pkin(1) * t92 + pkin(7)) * qJD(1);
t96 = sin(qJ(3));
t98 = cos(qJ(3));
t78 = t96 * qJD(2) + t98 * t83;
t93 = cos(pkin(9));
t84 = (-pkin(1) * t93 - pkin(2)) * qJD(1);
t105 = qJD(1) * t98;
t104 = t96 * qJD(1);
t103 = qJ(5) * qJD(1);
t75 = qJD(3) * qJ(4) + t78;
t76 = -pkin(3) * t105 - qJ(4) * t104 + t84;
t77 = t98 * qJD(2) - t96 * t83;
t102 = qJD(4) - t77;
t72 = pkin(4) * t105 + qJD(5) - t76;
t73 = t103 * t98 - t75;
t101 = -t103 * t96 + t102;
t97 = cos(qJ(6));
t95 = sin(qJ(6));
t85 = qJD(6) + t104;
t82 = -t95 * qJD(3) - t105 * t97;
t81 = -t97 * qJD(3) + t105 * t95;
t74 = -qJD(3) * pkin(3) + t102;
t71 = qJD(3) * pkin(5) - t73;
t70 = qJD(3) * t106 + t101;
t69 = (-pkin(8) + t106) * qJD(3) + t101;
t68 = (pkin(5) * t96 + pkin(8) * t98) * qJD(1) + t72;
t67 = t95 * t68 + t97 * t69;
t66 = t97 * t68 - t95 * t69;
t1 = qJD(2) ^ 2 * t107 + m(4) * (t77 ^ 2 + t78 ^ 2 + t84 ^ 2) / 0.2e1 + m(7) * (t66 ^ 2 + t67 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t70 ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t74 ^ 2 + t75 ^ 2 + t76 ^ 2) / 0.2e1 + (t66 * mrSges(7,1) - t67 * mrSges(7,2) + Ifges(7,3) * t85 / 0.2e1) * t85 + (t71 * mrSges(7,2) - t66 * mrSges(7,3) + Ifges(7,5) * t85 + Ifges(7,1) * t82 / 0.2e1) * t82 + (-t71 * mrSges(7,1) + t67 * mrSges(7,3) + Ifges(7,4) * t82 + Ifges(7,6) * t85 + Ifges(7,2) * t81 / 0.2e1) * t81 + (t77 * mrSges(4,1) - t74 * mrSges(5,1) - t73 * mrSges(6,1) - t78 * mrSges(4,2) + t70 * mrSges(6,2) + t75 * mrSges(5,3) + (Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * qJD(3)) * qJD(3) + ((-t84 * mrSges(4,1) - t76 * mrSges(5,1) + t75 * mrSges(5,2) - t72 * mrSges(6,2) + t78 * mrSges(4,3) + t73 * mrSges(6,3) + (Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1 + Ifges(6,1) / 0.2e1) * t105) * t98 + (-t70 * mrSges(6,3) - t77 * mrSges(4,3) + t74 * mrSges(5,2) + t84 * mrSges(4,2) + t72 * mrSges(6,1) - t76 * mrSges(5,3) + (Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1 + Ifges(6,2) / 0.2e1) * t104 + (Ifges(4,4) + Ifges(6,4) - Ifges(5,5)) * t105) * t96 + ((Ifges(6,5) + Ifges(4,6) - Ifges(5,6)) * t98 + (Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t96) * qJD(3) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t93 * mrSges(3,1) - t92 * mrSges(3,2) + (t92 ^ 2 + t93 ^ 2) * t107 * pkin(1)) * pkin(1)) * qJD(1)) * qJD(1);
T  = t1;
