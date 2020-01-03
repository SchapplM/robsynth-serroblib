% Calculate kinetic energy for
% S5RRRRR7
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d4,d5]';
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
% Datum: 2019-12-31 22:23
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR7_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR7_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 22:21:19
% EndTime: 2019-12-31 22:21:19
% DurationCPUTime: 0.45s
% Computational Cost: add. (474->86), mult. (1068->137), div. (0->0), fcn. (752->8), ass. (0->35)
t104 = qJD(1) * (-pkin(7) - pkin(6));
t102 = pkin(6) * mrSges(3,3);
t95 = sin(qJ(2));
t87 = qJD(2) * pkin(2) + t95 * t104;
t99 = cos(qJ(2));
t88 = t99 * t104;
t94 = sin(qJ(3));
t98 = cos(qJ(3));
t79 = t98 * t87 + t88 * t94;
t85 = (t94 * t99 + t95 * t98) * qJD(1);
t91 = qJD(2) + qJD(3);
t71 = pkin(3) * t91 - pkin(8) * t85 + t79;
t80 = t94 * t87 - t98 * t88;
t84 = (-t94 * t95 + t98 * t99) * qJD(1);
t73 = pkin(8) * t84 + t80;
t93 = sin(qJ(4));
t97 = cos(qJ(4));
t68 = t93 * t71 + t97 * t73;
t89 = (-pkin(2) * t99 - pkin(1)) * qJD(1);
t67 = t71 * t97 - t73 * t93;
t77 = t84 * t97 - t85 * t93;
t81 = -pkin(3) * t84 + t89;
t96 = cos(qJ(5));
t92 = sin(qJ(5));
t90 = qJD(4) + t91;
t78 = t84 * t93 + t85 * t97;
t76 = qJD(5) - t77;
t75 = t78 * t96 + t90 * t92;
t74 = -t78 * t92 + t90 * t96;
t69 = -pkin(4) * t77 - pkin(9) * t78 + t81;
t66 = pkin(9) * t90 + t68;
t65 = -pkin(4) * t90 - t67;
t64 = t66 * t96 + t69 * t92;
t63 = -t66 * t92 + t69 * t96;
t1 = Ifges(3,3) * qJD(2) ^ 2 / 0.2e1 + m(4) * (t79 ^ 2 + t80 ^ 2 + t89 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t81 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t65 ^ 2) / 0.2e1 + (t79 * mrSges(4,1) - t80 * mrSges(4,2) + Ifges(4,3) * t91 / 0.2e1) * t91 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * t90 / 0.2e1) * t90 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t76 / 0.2e1) * t76 + (t89 * mrSges(4,2) - t79 * mrSges(4,3) + Ifges(4,5) * t91 + Ifges(4,1) * t85 / 0.2e1) * t85 + (t81 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * t90 + Ifges(5,1) * t78 / 0.2e1) * t78 + (t65 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t76 + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t89 * mrSges(4,1) + t80 * mrSges(4,3) + Ifges(4,4) * t85 + Ifges(4,6) * t91 + Ifges(4,2) * t84 / 0.2e1) * t84 + (-t81 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t78 + Ifges(5,6) * t90 + Ifges(5,2) * t77 / 0.2e1) * t77 + (-t65 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,6) * t76 + Ifges(6,2) * t74 / 0.2e1) * t74 + ((m(3) * (pkin(1) ^ 2 + (t95 ^ 2 + t99 ^ 2) * pkin(6) ^ 2) / 0.2e1 + Ifges(2,3) / 0.2e1 + (pkin(1) * mrSges(3,1) + (t102 + Ifges(3,2) / 0.2e1) * t99) * t99 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t99 + (t102 + Ifges(3,1) / 0.2e1) * t95) * t95) * qJD(1) + ((-pkin(6) * mrSges(3,2) + Ifges(3,6)) * t99 + (-pkin(6) * mrSges(3,1) + Ifges(3,5)) * t95) * qJD(2)) * qJD(1);
T = t1;
