% Calculate kinetic energy for
% S5RPRRR9
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
% Datum: 2019-12-31 19:08
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR9_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR9_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:07:15
% EndTime: 2019-12-31 19:07:16
% DurationCPUTime: 0.44s
% Computational Cost: add. (413->80), mult. (1032->129), div. (0->0), fcn. (748->8), ass. (0->36)
t94 = cos(pkin(9));
t107 = t94 ^ 2;
t106 = qJD(1) * (pkin(6) + qJ(2));
t105 = m(3) / 0.2e1;
t100 = cos(qJ(3));
t93 = sin(pkin(9));
t86 = t93 * t106;
t87 = t94 * t106;
t97 = sin(qJ(3));
t78 = -t100 * t86 - t87 * t97;
t85 = (t100 * t93 + t94 * t97) * qJD(1);
t71 = qJD(3) * pkin(3) - pkin(7) * t85 + t78;
t79 = t100 * t87 - t97 * t86;
t84 = (t100 * t94 - t93 * t97) * qJD(1);
t72 = pkin(7) * t84 + t79;
t96 = sin(qJ(4));
t99 = cos(qJ(4));
t67 = t96 * t71 + t99 * t72;
t66 = t71 * t99 - t72 * t96;
t76 = t84 * t99 - t85 * t96;
t88 = qJD(2) + (-pkin(2) * t94 - pkin(1)) * qJD(1);
t80 = -pkin(3) * t84 + t88;
t98 = cos(qJ(5));
t95 = sin(qJ(5));
t92 = qJD(3) + qJD(4);
t89 = -qJD(1) * pkin(1) + qJD(2);
t77 = t84 * t96 + t85 * t99;
t75 = qJD(5) - t76;
t74 = t77 * t98 + t92 * t95;
t73 = -t77 * t95 + t92 * t98;
t68 = -pkin(4) * t76 - pkin(8) * t77 + t80;
t65 = pkin(8) * t92 + t67;
t64 = -pkin(4) * t92 - t66;
t63 = t65 * t98 + t68 * t95;
t62 = -t65 * t95 + t68 * t98;
t1 = t89 ^ 2 * t105 + m(4) * (t78 ^ 2 + t79 ^ 2 + t88 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t80 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * t92 / 0.2e1) * t92 + (t88 * mrSges(4,2) - t78 * mrSges(4,3) + Ifges(4,1) * t85 / 0.2e1) * t85 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t75 / 0.2e1) * t75 + (-t88 * mrSges(4,1) + t79 * mrSges(4,3) + Ifges(4,4) * t85 + Ifges(4,2) * t84 / 0.2e1) * t84 + (t80 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,5) * t92 + Ifges(5,1) * t77 / 0.2e1) * t77 + (t64 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t75 + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t80 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,6) * t92 + Ifges(5,2) * t76 / 0.2e1) * t76 + (-t64 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,6) * t75 + Ifges(6,2) * t73 / 0.2e1) * t73 + (t78 * mrSges(4,1) - t79 * mrSges(4,2) + Ifges(4,5) * t85 + Ifges(4,6) * t84 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t89 * (-mrSges(3,1) * t94 + mrSges(3,2) * t93) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t105 + mrSges(3,3)) * (t93 ^ 2 + t107) * qJ(2) + Ifges(3,2) * t107 / 0.2e1 + (Ifges(3,4) * t94 + Ifges(3,1) * t93 / 0.2e1) * t93) * qJD(1)) * qJD(1);
T = t1;
