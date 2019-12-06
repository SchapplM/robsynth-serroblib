% Calculate kinetic energy for
% S5RPRRR5
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
% Datum: 2019-12-05 18:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR5_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR5_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:16:04
% EndTime: 2019-12-05 18:16:04
% DurationCPUTime: 0.20s
% Computational Cost: add. (217->61), mult. (369->100), div. (0->0), fcn. (186->8), ass. (0->30)
t96 = m(3) / 0.2e1;
t82 = qJD(1) + qJD(3);
t89 = cos(qJ(4));
t95 = t82 * t89;
t84 = cos(pkin(9));
t76 = (pkin(1) * t84 + pkin(2)) * qJD(1);
t87 = sin(qJ(3));
t90 = cos(qJ(3));
t83 = sin(pkin(9));
t94 = pkin(1) * qJD(1) * t83;
t72 = t87 * t76 + t90 * t94;
t70 = t82 * pkin(7) + t72;
t86 = sin(qJ(4));
t66 = t86 * qJD(2) + t89 * t70;
t71 = t90 * t76 - t87 * t94;
t91 = qJD(2) ^ 2;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t81 = qJD(4) + qJD(5);
t80 = t89 * qJD(2);
t74 = (t85 * t89 + t86 * t88) * t82;
t73 = (-t85 * t86 + t88 * t89) * t82;
t69 = -t82 * pkin(3) - t71;
t67 = (-pkin(4) * t89 - pkin(3)) * t82 - t71;
t65 = -t86 * t70 + t80;
t64 = pkin(8) * t95 + t66;
t63 = qJD(4) * pkin(4) + t80 + (-pkin(8) * t82 - t70) * t86;
t62 = t85 * t63 + t88 * t64;
t61 = t88 * t63 - t85 * t64;
t1 = m(6) * (t61 ^ 2 + t62 ^ 2 + t67 ^ 2) / 0.2e1 + m(4) * (t71 ^ 2 + t72 ^ 2 + t91) / 0.2e1 + t91 * t96 + m(5) * (t65 ^ 2 + t66 ^ 2 + t69 ^ 2) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t65 * mrSges(5,1) - t66 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t67 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t74 / 0.2e1) * t74 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t84 * mrSges(3,1) - t83 * mrSges(3,2) + (t83 ^ 2 + t84 ^ 2) * t96 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (-t67 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,6) * t81 + Ifges(6,2) * t73 / 0.2e1) * t73 + (-t72 * mrSges(4,2) + t71 * mrSges(4,1) + Ifges(4,3) * t82 / 0.2e1 + (-t69 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t95 / 0.2e1) * t89 + (t69 * mrSges(5,2) - t65 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t89 + Ifges(5,1) * t86 / 0.2e1) * t82) * t86) * t82;
T = t1;
