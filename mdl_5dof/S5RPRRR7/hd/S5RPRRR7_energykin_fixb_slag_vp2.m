% Calculate kinetic energy for
% S5RPRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:04
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:03:08
% EndTime: 2019-12-31 19:03:08
% DurationCPUTime: 0.37s
% Computational Cost: add. (264->73), mult. (556->121), div. (0->0), fcn. (324->8), ass. (0->32)
t101 = m(3) / 0.2e1;
t89 = sin(pkin(9));
t84 = (pkin(1) * t89 + pkin(6)) * qJD(1);
t93 = sin(qJ(3));
t96 = cos(qJ(3));
t79 = t93 * qJD(2) + t96 * t84;
t76 = qJD(3) * pkin(7) + t79;
t90 = cos(pkin(9));
t99 = -pkin(1) * t90 - pkin(2);
t77 = (-pkin(3) * t96 - pkin(7) * t93 + t99) * qJD(1);
t92 = sin(qJ(4));
t95 = cos(qJ(4));
t68 = t95 * t76 + t92 * t77;
t100 = qJD(1) * t93;
t67 = -t76 * t92 + t95 * t77;
t78 = qJD(2) * t96 - t93 * t84;
t87 = -qJD(1) * t96 + qJD(4);
t75 = -qJD(3) * pkin(3) - t78;
t94 = cos(qJ(5));
t91 = sin(qJ(5));
t86 = qJD(5) + t87;
t85 = t99 * qJD(1);
t83 = qJD(3) * t92 + t95 * t100;
t82 = qJD(3) * t95 - t92 * t100;
t71 = t82 * t91 + t83 * t94;
t70 = t82 * t94 - t83 * t91;
t69 = -pkin(4) * t82 + t75;
t66 = pkin(8) * t82 + t68;
t65 = pkin(4) * t87 - pkin(8) * t83 + t67;
t64 = t65 * t91 + t66 * t94;
t63 = t65 * t94 - t66 * t91;
t1 = qJD(2) ^ 2 * t101 + m(4) * (t78 ^ 2 + t79 ^ 2 + t85 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t75 ^ 2) / 0.2e1 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t87 / 0.2e1) * t87 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t86 / 0.2e1) * t86 + (t75 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t87 + Ifges(5,1) * t83 / 0.2e1) * t83 + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t86 + Ifges(6,1) * t71 / 0.2e1) * t71 + (-t75 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t83 + Ifges(5,6) * t87 + Ifges(5,2) * t82 / 0.2e1) * t82 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t71 + Ifges(6,6) * t86 + Ifges(6,2) * t70 / 0.2e1) * t70 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t85 * (-mrSges(4,1) * t96 + mrSges(4,2) * t93) + (-t78 * t93 + t79 * t96) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t93 + Ifges(4,6) * t96) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t90 * mrSges(3,1) - t89 * mrSges(3,2) + (t89 ^ 2 + t90 ^ 2) * t101 * pkin(1)) * pkin(1) + Ifges(4,2) * t96 ^ 2 / 0.2e1 + (Ifges(4,4) * t96 + Ifges(4,1) * t93 / 0.2e1) * t93) * qJD(1)) * qJD(1);
T = t1;
