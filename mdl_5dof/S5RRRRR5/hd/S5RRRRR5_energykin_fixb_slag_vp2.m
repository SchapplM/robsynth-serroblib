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
% Datum: 2019-12-05 18:59
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:58:07
% EndTime: 2019-12-05 18:58:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (291->56), mult. (357->99), div. (0->0), fcn. (174->8), ass. (0->29)
t80 = qJD(1) + qJD(2);
t78 = qJD(3) + t80;
t95 = t78 / 0.2e1;
t88 = cos(qJ(2));
t93 = qJD(1) * pkin(1);
t75 = t80 * pkin(2) + t88 * t93;
t83 = sin(qJ(3));
t87 = cos(qJ(3));
t84 = sin(qJ(2));
t92 = t84 * t93;
t71 = t83 * t75 + t87 * t92;
t69 = t78 * pkin(8) + t71;
t94 = t69 * mrSges(5,3);
t91 = pkin(9) * t78 + t69;
t70 = t87 * t75 - t83 * t92;
t86 = cos(qJ(4));
t85 = cos(qJ(5));
t82 = sin(qJ(4));
t81 = sin(qJ(5));
t79 = qJD(4) + qJD(5);
t73 = (t81 * t86 + t82 * t85) * t78;
t72 = (-t81 * t82 + t85 * t86) * t78;
t68 = -t78 * pkin(3) - t70;
t66 = (-pkin(4) * t86 - pkin(3)) * t78 - t70;
t65 = t91 * t86;
t64 = qJD(4) * pkin(4) - t91 * t82;
t63 = t81 * t64 + t85 * t65;
t62 = t85 * t64 - t81 * t65;
t1 = m(6) * (t62 ^ 2 + t63 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t68 ^ 2 + (t82 ^ 2 + t86 ^ 2) * t69 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t84 ^ 2 + t88 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (Ifges(3,3) * t80 / 0.2e1 + (mrSges(3,1) * t88 - mrSges(3,2) * t84) * t93) * t80 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t79 / 0.2e1) * t79 + (Ifges(5,3) * qJD(4) / 0.2e1 + (-t82 * mrSges(5,1) - t86 * mrSges(5,2)) * t69) * qJD(4) + (t66 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t79 + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t66 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,6) * t79 + Ifges(6,2) * t72 / 0.2e1) * t72 + (-t71 * mrSges(4,2) + t70 * mrSges(4,1) + Ifges(4,3) * t95 + (-t68 * mrSges(5,1) + Ifges(5,6) * qJD(4) + (Ifges(5,2) * t95 + t94) * t86) * t86 + (Ifges(5,4) * t86 * t78 + t68 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t95 + t94) * t82) * t82) * t78;
T = t1;
