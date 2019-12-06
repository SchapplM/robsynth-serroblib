% Calculate kinetic energy for
% S5RRPRR4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d4,d5,theta3]';
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
% Datum: 2019-12-05 18:32
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPRR4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPRR4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:31:47
% EndTime: 2019-12-05 18:31:48
% DurationCPUTime: 0.20s
% Computational Cost: add. (239->60), mult. (368->101), div. (0->0), fcn. (186->8), ass. (0->29)
t83 = qJD(1) + qJD(2);
t90 = cos(qJ(4));
t96 = t83 * t90;
t91 = cos(qJ(2));
t95 = qJD(1) * pkin(1);
t76 = pkin(2) * t83 + t91 * t95;
t84 = sin(pkin(9));
t85 = cos(pkin(9));
t88 = sin(qJ(2));
t94 = t88 * t95;
t72 = t84 * t76 + t85 * t94;
t70 = pkin(7) * t83 + t72;
t87 = sin(qJ(4));
t66 = t87 * qJD(3) + t90 * t70;
t71 = t85 * t76 - t84 * t94;
t89 = cos(qJ(5));
t86 = sin(qJ(5));
t82 = qJD(4) + qJD(5);
t80 = t90 * qJD(3);
t74 = (t86 * t90 + t87 * t89) * t83;
t73 = (-t86 * t87 + t89 * t90) * t83;
t69 = -pkin(3) * t83 - t71;
t67 = (-pkin(4) * t90 - pkin(3)) * t83 - t71;
t65 = -t70 * t87 + t80;
t64 = pkin(8) * t96 + t66;
t63 = qJD(4) * pkin(4) + t80 + (-pkin(8) * t83 - t70) * t87;
t62 = t63 * t86 + t64 * t89;
t61 = t63 * t89 - t64 * t86;
t1 = m(6) * (t61 ^ 2 + t62 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t71 ^ 2 + t72 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t88 ^ 2 + t91 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t65 * mrSges(5,1) - t66 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t67 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t67 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,6) * t82 + Ifges(6,2) * t73 / 0.2e1) * t73 + (t71 * mrSges(4,1) - t72 * mrSges(4,2) + (mrSges(3,1) * t91 - mrSges(3,2) * t88) * t95 + (-t69 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t96 / 0.2e1) * t90 + (t69 * mrSges(5,2) - t65 * mrSges(5,3) + Ifges(5,5) * qJD(4)) * t87 + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t90 + Ifges(5,1) * t87 / 0.2e1) * t87) * t83) * t83;
T = t1;
