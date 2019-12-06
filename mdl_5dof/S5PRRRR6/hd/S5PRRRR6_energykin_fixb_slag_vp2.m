% Calculate kinetic energy for
% S5PRRRR6
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d2,d3,d4,d5,theta1]';
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
% Datum: 2019-12-05 17:10
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PRRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PRRRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:09:30
% EndTime: 2019-12-05 17:09:30
% DurationCPUTime: 0.18s
% Computational Cost: add. (203->55), mult. (311->96), div. (0->0), fcn. (174->8), ass. (0->27)
t81 = qJD(2) + qJD(3);
t94 = t81 / 0.2e1;
t89 = cos(qJ(2));
t77 = qJD(2) * pkin(2) + t89 * qJD(1);
t84 = sin(qJ(3));
t88 = cos(qJ(3));
t85 = sin(qJ(2));
t92 = qJD(1) * t85;
t73 = t84 * t77 + t88 * t92;
t71 = t81 * pkin(7) + t73;
t93 = t71 * mrSges(5,3);
t91 = pkin(8) * t81 + t71;
t72 = t88 * t77 - t84 * t92;
t87 = cos(qJ(4));
t86 = cos(qJ(5));
t83 = sin(qJ(4));
t82 = sin(qJ(5));
t80 = qJD(4) + qJD(5);
t75 = (t82 * t87 + t83 * t86) * t81;
t74 = (-t82 * t83 + t86 * t87) * t81;
t70 = -t81 * pkin(3) - t72;
t68 = (-pkin(4) * t87 - pkin(3)) * t81 - t72;
t67 = t91 * t87;
t66 = qJD(4) * pkin(4) - t91 * t83;
t65 = t82 * t66 + t86 * t67;
t64 = t86 * t66 - t82 * t67;
t1 = m(6) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t70 ^ 2 + (t83 ^ 2 + t87 ^ 2) * t71 ^ 2) / 0.2e1 + m(4) * (t72 ^ 2 + t73 ^ 2) / 0.2e1 + (m(3) * (t85 ^ 2 + t89 ^ 2) / 0.2e1 + m(2) / 0.2e1) * qJD(1) ^ 2 + (t64 * mrSges(6,1) - t65 * mrSges(6,2) + Ifges(6,3) * t80 / 0.2e1) * t80 + (Ifges(5,3) * qJD(4) / 0.2e1 + (-t83 * mrSges(5,1) - t87 * mrSges(5,2)) * t71) * qJD(4) + (Ifges(3,3) * qJD(2) / 0.2e1 + (t89 * mrSges(3,1) - t85 * mrSges(3,2)) * qJD(1)) * qJD(2) + (t68 * mrSges(6,2) - t64 * mrSges(6,3) + Ifges(6,5) * t80 + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t68 * mrSges(6,1) + t65 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,6) * t80 + Ifges(6,2) * t74 / 0.2e1) * t74 + (-t73 * mrSges(4,2) + t72 * mrSges(4,1) + Ifges(4,3) * t94 + (-t70 * mrSges(5,1) + Ifges(5,6) * qJD(4) + (Ifges(5,2) * t94 + t93) * t87) * t87 + (Ifges(5,4) * t87 * t81 + t70 * mrSges(5,2) + Ifges(5,5) * qJD(4) + (Ifges(5,1) * t94 + t93) * t83) * t83) * t81;
T = t1;
