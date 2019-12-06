% Calculate kinetic energy for
% S5RPPPR1
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d5,theta2,theta3,theta4]';
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
% Datum: 2019-12-05 17:29
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPPPR1_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPPPR1_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:28:43
% EndTime: 2019-12-05 17:28:43
% DurationCPUTime: 0.32s
% Computational Cost: add. (197->66), mult. (483->110), div. (0->0), fcn. (276->8), ass. (0->30)
t94 = m(3) / 0.2e1;
t81 = sin(pkin(8));
t84 = cos(pkin(8));
t85 = cos(pkin(7));
t90 = -pkin(1) * t85 - pkin(2);
t69 = qJD(3) + (-pkin(3) * t84 - qJ(4) * t81 + t90) * qJD(1);
t82 = sin(pkin(7));
t78 = (pkin(1) * t82 + qJ(3)) * qJD(1);
t75 = qJD(2) * t81 + t78 * t84;
t80 = sin(pkin(9));
t83 = cos(pkin(9));
t65 = t80 * t69 + t83 * t75;
t93 = qJD(1) * t81;
t92 = qJD(1) * t84;
t91 = t80 * t93;
t64 = t83 * t69 - t75 * t80;
t74 = qJD(2) * t84 - t81 * t78;
t73 = qJD(4) - t74;
t87 = cos(qJ(5));
t86 = sin(qJ(5));
t79 = qJD(5) - t92;
t77 = t90 * qJD(1) + qJD(3);
t72 = (-t80 * t86 + t83 * t87) * t93;
t71 = (-t80 * t87 - t83 * t86) * t93;
t66 = pkin(4) * t91 + t73;
t63 = -pkin(6) * t91 + t65;
t62 = (-pkin(6) * t81 * t83 - pkin(4) * t84) * qJD(1) + t64;
t61 = t62 * t86 + t63 * t87;
t60 = t62 * t87 - t63 * t86;
t1 = qJD(2) ^ 2 * t94 + m(6) * (t60 ^ 2 + t61 ^ 2 + t66 ^ 2) / 0.2e1 + m(5) * (t64 ^ 2 + t65 ^ 2 + t73 ^ 2) / 0.2e1 + m(4) * (t74 ^ 2 + t75 ^ 2 + t77 ^ 2) / 0.2e1 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,3) * t79 / 0.2e1) * t79 + (t66 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,5) * t79 + Ifges(6,1) * t72 / 0.2e1) * t72 + (-t66 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t72 + Ifges(6,6) * t79 + Ifges(6,2) * t71 / 0.2e1) * t71 + ((t73 * (mrSges(5,1) * t80 + mrSges(5,2) * t83) + t77 * mrSges(4,2) - t74 * mrSges(4,3) + (Ifges(5,1) * t83 ^ 2 / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t83 + Ifges(5,2) * t80 / 0.2e1) * t80) * t93 + (-t64 * t83 - t65 * t80) * mrSges(5,3)) * t81 + (t65 * mrSges(5,2) - t77 * mrSges(4,1) - t64 * mrSges(5,1) + t75 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t92 + (-Ifges(5,5) * t83 + Ifges(5,6) * t80 + Ifges(4,4)) * t93) * t84 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t85 * mrSges(3,1) - t82 * mrSges(3,2) + (t82 ^ 2 + t85 ^ 2) * t94 * pkin(1)) * pkin(1)) * qJD(1)) * qJD(1);
T = t1;
