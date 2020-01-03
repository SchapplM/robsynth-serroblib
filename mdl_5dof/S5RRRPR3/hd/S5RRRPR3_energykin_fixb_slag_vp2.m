% Calculate kinetic energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2020-01-03 12:10
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:08:47
% EndTime: 2020-01-03 12:08:48
% DurationCPUTime: 0.25s
% Computational Cost: add. (400->71), mult. (579->119), div. (0->0), fcn. (346->8), ass. (0->32)
t84 = qJD(1) + qJD(2);
t99 = t84 / 0.2e1;
t89 = sin(qJ(2));
t97 = qJD(1) * pkin(1);
t81 = pkin(7) * t84 + t89 * t97;
t98 = t81 * mrSges(4,3);
t88 = sin(qJ(3));
t95 = qJ(4) * t84 + t81;
t75 = qJD(3) * pkin(3) - t95 * t88;
t91 = cos(qJ(3));
t76 = t95 * t91;
t85 = sin(pkin(9));
t86 = cos(pkin(9));
t68 = t85 * t75 + t86 * t76;
t92 = cos(qJ(2));
t96 = t92 * t97;
t67 = t86 * t75 - t76 * t85;
t77 = -t96 + qJD(4) + (-pkin(3) * t91 - pkin(2)) * t84;
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t83 = qJD(3) + qJD(5);
t82 = -pkin(2) * t84 - t96;
t79 = (t85 * t91 + t86 * t88) * t84;
t78 = (-t85 * t88 + t86 * t91) * t84;
t71 = -pkin(4) * t78 + t77;
t70 = t78 * t87 + t79 * t90;
t69 = t78 * t90 - t79 * t87;
t66 = pkin(8) * t78 + t68;
t65 = qJD(3) * pkin(4) - pkin(8) * t79 + t67;
t64 = t65 * t87 + t66 * t90;
t63 = t65 * t90 - t66 * t87;
t1 = m(4) * (t82 ^ 2 + (t88 ^ 2 + t91 ^ 2) * t81 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t77 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t89 ^ 2 + t92 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (t77 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,1) * t79 / 0.2e1) * t79 + (-t77 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,2) * t78 / 0.2e1) * t78 + (t71 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t71 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t83 + Ifges(6,2) * t69 / 0.2e1) * t69 + (Ifges(3,3) * t99 + (mrSges(3,1) * t92 - mrSges(3,2) * t89) * t97 + (-t82 * mrSges(4,1) + (Ifges(4,2) * t99 + t98) * t91) * t91 + (Ifges(4,4) * t91 * t84 + t82 * mrSges(4,2) + (Ifges(4,1) * t99 + t98) * t88) * t88) * t84 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,5) * t79 + Ifges(5,6) * t78 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t88 + Ifges(4,6) * t91) * t84 + (-mrSges(4,1) * t88 - mrSges(4,2) * t91) * t81) * qJD(3);
T = t1;
