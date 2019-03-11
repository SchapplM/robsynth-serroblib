% Calculate kinetic energy for
% S6RPRPPR8
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% qJD [6x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d6]';
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
% Datum: 2019-03-09 03:00
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S6RPRPPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(6,1),zeros(8,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: qJ has to be [6x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [6 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: qJD has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRPPR8_energykin_fixb_slag_vp2: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-03-09 02:58:41
% EndTime: 2019-03-09 02:58:42
% DurationCPUTime: 0.27s
% Computational Cost: add. (217->92), mult. (394->120), div. (0->0), fcn. (144->4), ass. (0->33)
t94 = -pkin(3) - pkin(4);
t75 = qJD(2) + (-pkin(1) - pkin(7)) * qJD(1);
t93 = t75 * mrSges(4,3);
t82 = sin(qJ(3));
t70 = qJD(3) * qJ(4) + t75 * t82;
t92 = t82 * qJD(1);
t84 = cos(qJ(3));
t91 = t84 * qJD(1);
t77 = qJ(4) * t91;
t90 = qJD(5) + t77;
t89 = qJ(5) * qJD(1);
t88 = qJD(1) * qJ(2);
t87 = -pkin(8) + t94;
t67 = -t82 * t89 - t70;
t86 = qJD(4) + (-t75 - t89) * t84;
t85 = qJD(1) ^ 2;
t83 = cos(qJ(6));
t81 = sin(qJ(6));
t80 = t85 * qJ(2) ^ 2;
t78 = -pkin(1) * qJD(1) + qJD(2);
t76 = qJD(6) + t91;
t72 = -qJD(3) * t81 + t83 * t92;
t71 = -qJD(3) * t83 - t81 * t92;
t69 = -t77 + (pkin(3) * t82 + qJ(2)) * qJD(1);
t68 = -qJD(3) * pkin(3) - t75 * t84 + qJD(4);
t66 = (t82 * t94 - qJ(2)) * qJD(1) + t90;
t65 = qJD(3) * pkin(5) - t67;
t64 = qJD(3) * t94 + t86;
t63 = qJD(3) * t87 + t86;
t62 = (pkin(5) * t84 + t82 * t87 - qJ(2)) * qJD(1) + t90;
t61 = t62 * t81 + t63 * t83;
t60 = t62 * t83 - t63 * t81;
t1 = m(4) * (t80 + (t82 ^ 2 + t84 ^ 2) * t75 ^ 2) / 0.2e1 + m(3) * (t78 ^ 2 + t80) / 0.2e1 + m(7) * (t60 ^ 2 + t61 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t66 ^ 2 + t67 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + t69 ^ 2 + t70 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + Ifges(3,1) / 0.2e1 + qJ(2) * mrSges(3,3)) * t85 + (t60 * mrSges(7,1) - t61 * mrSges(7,2) + Ifges(7,3) * t76 / 0.2e1) * t76 + (t65 * mrSges(7,2) - t60 * mrSges(7,3) + Ifges(7,5) * t76 + Ifges(7,1) * t72 / 0.2e1) * t72 + (-t65 * mrSges(7,1) + t61 * mrSges(7,3) + Ifges(7,4) * t72 + Ifges(7,6) * t76 + Ifges(7,2) * t71 / 0.2e1) * t71 + (t78 * mrSges(3,2) + (mrSges(4,2) * t88 + t66 * mrSges(6,1) + t68 * mrSges(5,2) - t69 * mrSges(5,3) - t64 * mrSges(6,3) + (-t93 + (Ifges(6,2) / 0.2e1 + Ifges(5,1) / 0.2e1 + Ifges(4,1) / 0.2e1) * qJD(1)) * t84) * t84 + (t66 * mrSges(6,2) + t69 * mrSges(5,1) - t70 * mrSges(5,2) - t67 * mrSges(6,3) + mrSges(4,1) * t88 + (-t93 + (Ifges(6,1) / 0.2e1 + Ifges(5,3) / 0.2e1 + Ifges(4,2) / 0.2e1) * qJD(1)) * t82 + (-Ifges(4,4) - Ifges(6,4) + Ifges(5,5)) * t91) * t82) * qJD(1) + (-t68 * mrSges(5,1) - t67 * mrSges(6,1) + t64 * mrSges(6,2) + t70 * mrSges(5,3) + (Ifges(6,3) / 0.2e1 + Ifges(5,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * qJD(3) + (mrSges(4,1) * t84 - mrSges(4,2) * t82) * t75 + ((Ifges(5,4) + Ifges(4,5) + Ifges(6,6)) * t84 + (-Ifges(6,5) - Ifges(4,6) + Ifges(5,6)) * t82) * qJD(1)) * qJD(3);
T  = t1;
