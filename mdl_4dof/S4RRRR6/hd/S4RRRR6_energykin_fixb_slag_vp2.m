% Calculate kinetic energy for
% S4RRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [4x1]
%   Generalized joint coordinates (joint angles)
% qJD [4x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,alpha2,d1,d2,d3,d4]';
% m_mdh [5x1]
%   mass of all robot links (including the base)
% mrSges [5x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [5x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% T [1x1]
%   kinetic energy

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-31 17:31
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S4RRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(4,1),zeros(4,1),zeros(8,1),zeros(5,1),zeros(5,3),zeros(5,6)}
assert(isreal(qJ) && all(size(qJ) == [4 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: qJ has to be [4x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [4 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: qJD has to be [4x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [5 1]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: m has to be [5x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [5,3]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: mrSges has to be [5x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [5 6]), ...
  'S4RRRR6_energykin_fixb_slag_vp2: Ifges has to be [5x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 17:29:19
% EndTime: 2019-12-31 17:29:19
% DurationCPUTime: 0.25s
% Computational Cost: add. (278->66), mult. (662->113), div. (0->0), fcn. (468->8), ass. (0->33)
t89 = sin(qJ(2));
t92 = cos(qJ(2));
t85 = sin(pkin(4));
t98 = qJD(1) * t85;
t94 = t92 * t98;
t97 = cos(pkin(4)) * qJD(1);
t96 = pkin(1) * t97;
t79 = pkin(6) * t94 + t89 * t96;
t84 = qJD(2) + t97;
t73 = pkin(7) * t84 + t79;
t74 = (-pkin(2) * t92 - pkin(7) * t89 - pkin(1)) * t98;
t88 = sin(qJ(3));
t91 = cos(qJ(3));
t67 = t91 * t73 + t88 * t74;
t95 = t89 * t98;
t66 = -t73 * t88 + t74 * t91;
t78 = -pkin(6) * t95 + t92 * t96;
t76 = t84 * t91 - t88 * t95;
t72 = -pkin(2) * t84 - t78;
t93 = qJD(1) ^ 2;
t90 = cos(qJ(4));
t87 = sin(qJ(4));
t80 = qJD(3) - t94;
t77 = t84 * t88 + t91 * t95;
t75 = qJD(4) - t76;
t69 = t77 * t90 + t80 * t87;
t68 = -t77 * t87 + t80 * t90;
t65 = pkin(8) * t80 + t67;
t64 = -pkin(3) * t80 - t66;
t63 = -pkin(3) * t76 - pkin(8) * t77 + t72;
t62 = t63 * t87 + t65 * t90;
t61 = t63 * t90 - t65 * t87;
t1 = m(3) * (pkin(1) ^ 2 * t85 ^ 2 * t93 + t78 ^ 2 + t79 ^ 2) / 0.2e1 + t93 * Ifges(2,3) / 0.2e1 + m(5) * (t61 ^ 2 + t62 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t72 ^ 2) / 0.2e1 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + Ifges(4,3) * t80 / 0.2e1) * t80 + (t61 * mrSges(5,1) - t62 * mrSges(5,2) + Ifges(5,3) * t75 / 0.2e1) * t75 + (t72 * mrSges(4,2) - t66 * mrSges(4,3) + Ifges(4,5) * t80 + Ifges(4,1) * t77 / 0.2e1) * t77 + (t64 * mrSges(5,2) - t61 * mrSges(5,3) + Ifges(5,5) * t75 + Ifges(5,1) * t69 / 0.2e1) * t69 + (((pkin(1) * mrSges(3,1) + Ifges(3,2) * t92 / 0.2e1) * t92 + (-pkin(1) * mrSges(3,2) + Ifges(3,4) * t92 + Ifges(3,1) * t89 / 0.2e1) * t89) * t98 + (-t78 * t89 + t79 * t92) * mrSges(3,3)) * t98 + (t78 * mrSges(3,1) - t79 * mrSges(3,2) + Ifges(3,3) * t84 / 0.2e1 + (Ifges(3,5) * t89 + Ifges(3,6) * t92) * t98) * t84 + (-t72 * mrSges(4,1) + t67 * mrSges(4,3) + Ifges(4,4) * t77 + Ifges(4,6) * t80 + Ifges(4,2) * t76 / 0.2e1) * t76 + (-t64 * mrSges(5,1) + t62 * mrSges(5,3) + Ifges(5,4) * t69 + Ifges(5,6) * t75 + Ifges(5,2) * t68 / 0.2e1) * t68;
T = t1;
