% Calculate kinetic energy for
% S5RRPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d5,theta3,theta4]';
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
% Datum: 2019-12-05 18:19
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRPPR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:18:08
% EndTime: 2019-12-05 18:18:08
% DurationCPUTime: 0.16s
% Computational Cost: add. (224->54), mult. (353->93), div. (0->0), fcn. (184->8), ass. (0->29)
t81 = qJD(1) + qJD(2);
t95 = pkin(7) * t81;
t89 = cos(qJ(2));
t94 = qJD(1) * pkin(1);
t75 = t81 * pkin(2) + t89 * t94;
t83 = sin(pkin(8));
t85 = cos(pkin(8));
t87 = sin(qJ(2));
t93 = t87 * t94;
t71 = t83 * t75 + t85 * t93;
t69 = t81 * qJ(4) + t71;
t82 = sin(pkin(9));
t84 = cos(pkin(9));
t65 = t82 * qJD(3) + t84 * t69;
t70 = t85 * t75 - t83 * t93;
t92 = qJD(4) - t70;
t88 = cos(qJ(5));
t86 = sin(qJ(5));
t79 = t84 * qJD(3);
t73 = (t82 * t88 + t84 * t86) * t81;
t72 = (-t82 * t86 + t84 * t88) * t81;
t68 = -t81 * pkin(3) + t92;
t66 = (-pkin(4) * t84 - pkin(3)) * t81 + t92;
t64 = -t82 * t69 + t79;
t63 = t84 * t95 + t65;
t62 = t79 + (-t69 - t95) * t82;
t61 = t86 * t62 + t88 * t63;
t60 = t88 * t62 - t86 * t63;
t1 = m(6) * (t60 ^ 2 + t61 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (qJD(3) ^ 2 + t70 ^ 2 + t71 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t87 ^ 2 + t89 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t66 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t66 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,2) * t72 / 0.2e1) * t72 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,5) * t73 + Ifges(6,6) * t72 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t71 * mrSges(4,2) + t70 * mrSges(4,1) + t68 * (-mrSges(5,1) * t84 + mrSges(5,2) * t82) + (mrSges(3,1) * t89 - mrSges(3,2) * t87) * t94 + (-t64 * t82 + t65 * t84) * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) * t84 ^ 2 / 0.2e1 + (Ifges(5,4) * t84 + Ifges(5,1) * t82 / 0.2e1) * t82) * t81) * t81;
T = t1;
