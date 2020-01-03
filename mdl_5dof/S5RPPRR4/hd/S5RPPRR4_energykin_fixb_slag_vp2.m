% Calculate kinetic energy for
% S5RPPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d4,d5,theta2,theta3]';
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
% Datum: 2020-01-03 11:32
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:30:21
% EndTime: 2020-01-03 11:30:22
% DurationCPUTime: 0.47s
% Computational Cost: add. (335->80), mult. (898->132), div. (0->0), fcn. (604->8), ass. (0->37)
t107 = m(3) / 0.2e1;
t92 = sin(pkin(8));
t94 = cos(pkin(8));
t80 = qJD(2) + (-pkin(2) * t94 - qJ(3) * t92 - pkin(1)) * qJD(1);
t93 = cos(pkin(9));
t79 = t93 * t80;
t91 = sin(pkin(9));
t70 = t79 + (-pkin(6) * t92 * t93 + (-qJ(2) * t91 - pkin(3)) * t94) * qJD(1);
t105 = t92 * qJD(1);
t102 = t91 * t105;
t103 = qJD(1) * qJ(2);
t101 = t94 * t103;
t75 = t93 * t101 + t91 * t80;
t73 = -pkin(6) * t102 + t75;
t96 = sin(qJ(4));
t98 = cos(qJ(4));
t65 = t96 * t70 + t98 * t73;
t104 = t94 * qJD(1);
t85 = t92 * t103 + qJD(3);
t81 = pkin(3) * t102 + t85;
t64 = t98 * t70 - t96 * t73;
t86 = qJD(4) - t104;
t97 = cos(qJ(5));
t95 = sin(qJ(5));
t88 = -qJD(1) * pkin(1) + qJD(2);
t84 = qJD(5) + t86;
t77 = (-t91 * t96 + t93 * t98) * t105;
t76 = (-t91 * t98 - t93 * t96) * t105;
t74 = -t101 * t91 + t79;
t72 = -t76 * pkin(4) + t81;
t67 = t95 * t76 + t97 * t77;
t66 = t97 * t76 - t95 * t77;
t63 = t76 * pkin(7) + t65;
t62 = t86 * pkin(4) - t77 * pkin(7) + t64;
t61 = t95 * t62 + t97 * t63;
t60 = t97 * t62 - t95 * t63;
t1 = m(6) * (t60 ^ 2 + t61 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t81 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t85 ^ 2) / 0.2e1 + t88 ^ 2 * t107 + (t64 * mrSges(5,1) - t65 * mrSges(5,2) + Ifges(5,3) * t86 / 0.2e1) * t86 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t84 / 0.2e1) * t84 + (t81 * mrSges(5,2) - t64 * mrSges(5,3) + Ifges(5,5) * t86 + Ifges(5,1) * t77 / 0.2e1) * t77 + (t72 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t84 + Ifges(6,1) * t67 / 0.2e1) * t67 + (-t81 * mrSges(5,1) + t65 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,6) * t86 + Ifges(5,2) * t76 / 0.2e1) * t76 + (-t72 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t67 + Ifges(6,6) * t84 + Ifges(6,2) * t66 / 0.2e1) * t66 + ((t85 * (mrSges(4,1) * t91 + mrSges(4,2) * t93) + t88 * mrSges(3,2) + (Ifges(3,1) / 0.2e1 + Ifges(4,1) * t93 ^ 2 / 0.2e1 + (-Ifges(4,4) * t93 + Ifges(4,2) * t91 / 0.2e1) * t91) * t105 + (-t74 * t93 - t75 * t91) * mrSges(4,3)) * t92 + (-t88 * mrSges(3,1) + t75 * mrSges(4,2) - t74 * mrSges(4,1) + (Ifges(3,2) / 0.2e1 + Ifges(4,3) / 0.2e1) * t104 + (-Ifges(4,5) * t93 + Ifges(4,6) * t91 + Ifges(3,4)) * t105) * t94 + (Ifges(2,3) / 0.2e1 + (qJ(2) * t107 + mrSges(3,3)) * (t92 ^ 2 + t94 ^ 2) * qJ(2)) * qJD(1)) * qJD(1);
T = t1;
