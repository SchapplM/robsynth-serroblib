% Calculate kinetic energy for
% S5RRRRR5
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
% Datum: 2020-01-03 12:13
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRRR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 12:13:01
% EndTime: 2020-01-03 12:13:01
% DurationCPUTime: 0.16s
% Computational Cost: add. (291->56), mult. (357->99), div. (0->0), fcn. (174->8), ass. (0->29)
t82 = qJD(1) + qJD(2);
t80 = qJD(3) + t82;
t97 = t80 / 0.2e1;
t90 = cos(qJ(2));
t95 = qJD(1) * pkin(1);
t77 = t82 * pkin(2) + t90 * t95;
t85 = sin(qJ(3));
t89 = cos(qJ(3));
t86 = sin(qJ(2));
t94 = t86 * t95;
t73 = t85 * t77 + t89 * t94;
t71 = t80 * pkin(8) + t73;
t96 = t71 * mrSges(5,3);
t93 = pkin(9) * t80 + t71;
t72 = t89 * t77 - t85 * t94;
t88 = cos(qJ(4));
t87 = cos(qJ(5));
t84 = sin(qJ(4));
t83 = sin(qJ(5));
t81 = qJD(4) + qJD(5);
t75 = (t83 * t88 + t84 * t87) * t80;
t74 = (-t83 * t84 + t87 * t88) * t80;
t70 = -t80 * pkin(3) - t72;
t68 = (-pkin(4) * t88 - pkin(3)) * t80 - t72;
t67 = t93 * t88;
t66 = qJD(4) * pkin(4) - t93 * t84;
t65 = t83 * t66 + t87 * t67;
t64 = t87 * t66 - t83 * t67;
t1 = m(5) * (t70 ^ 2 + (t84 ^ 2 + t88 ^ 2) * t71 ^ 2) / 0.2e1 + m(6) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + m(4) * (t72 ^ 2 + t73 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t86 ^ 2 + t90 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t82 / 0.2e1 + (mrSges(3,1) * t90 - mrSges(3,2) * t86) * t95) * t82 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (Ifges(5,3) * qJD(4) / 0.2e1 + (-t84 * mrSges(5,1) - t88 * mrSges(5,2)) * t71) * qJD(4) + (t68 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t68 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,6) * t81 + Ifges(6,2) * t74 / 0.2e1) * t74 + (-t73 * mrSges(4,2) + t72 * mrSges(4,1) + Ifges(4,3) * t97 + (-t70 * mrSges(5,1) + Ifges(5,6) * qJD(4) + (Ifges(5,2) * t97 + t96) * t88) * t88 + (Ifges(5,4) * t88 * t80 + t70 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t97 + t96) * t84) * t84) * t80;
T = t1;
