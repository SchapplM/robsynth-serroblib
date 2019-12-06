% Calculate kinetic energy for
% S5RPRPR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2,theta4]';
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
% Datum: 2019-12-05 17:50
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 17:49:18
% EndTime: 2019-12-05 17:49:18
% DurationCPUTime: 0.17s
% Computational Cost: add. (202->55), mult. (354->92), div. (0->0), fcn. (184->8), ass. (0->30)
t95 = m(3) / 0.2e1;
t80 = qJD(1) + qJD(3);
t94 = pkin(7) * t80;
t84 = cos(pkin(8));
t75 = (pkin(1) * t84 + pkin(2)) * qJD(1);
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t82 = sin(pkin(8));
t93 = pkin(1) * qJD(1) * t82;
t71 = t86 * t75 + t88 * t93;
t69 = t80 * qJ(4) + t71;
t81 = sin(pkin(9));
t83 = cos(pkin(9));
t65 = t81 * qJD(2) + t83 * t69;
t70 = t88 * t75 - t86 * t93;
t92 = qJD(4) - t70;
t89 = qJD(2) ^ 2;
t87 = cos(qJ(5));
t85 = sin(qJ(5));
t79 = t83 * qJD(2);
t73 = (t81 * t87 + t83 * t85) * t80;
t72 = (-t81 * t85 + t83 * t87) * t80;
t68 = -t80 * pkin(3) + t92;
t66 = (-pkin(4) * t83 - pkin(3)) * t80 + t92;
t64 = -t81 * t69 + t79;
t63 = t83 * t94 + t65;
t62 = t79 + (-t69 - t94) * t81;
t61 = t85 * t62 + t87 * t63;
t60 = t87 * t62 - t85 * t63;
t1 = m(6) * (t60 ^ 2 + t61 ^ 2 + t66 ^ 2) / 0.2e1 + m(4) * (t70 ^ 2 + t71 ^ 2 + t89) / 0.2e1 + t89 * t95 + m(5) * (t64 ^ 2 + t65 ^ 2 + t68 ^ 2) / 0.2e1 + (t66 * mrSges(6,2) - t60 * mrSges(6,3) + Ifges(6,1) * t73 / 0.2e1) * t73 + (-t66 * mrSges(6,1) + t61 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,2) * t72 / 0.2e1) * t72 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t84 * mrSges(3,1) - t82 * mrSges(3,2) + (t82 ^ 2 + t84 ^ 2) * t95 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (t60 * mrSges(6,1) - t61 * mrSges(6,2) + Ifges(6,5) * t73 + Ifges(6,6) * t72 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (t68 * (-mrSges(5,1) * t83 + mrSges(5,2) * t81) - t71 * mrSges(4,2) + t70 * mrSges(4,1) + (Ifges(5,2) * t83 ^ 2 / 0.2e1 + Ifges(4,3) / 0.2e1 + (Ifges(5,4) * t83 + Ifges(5,1) * t81 / 0.2e1) * t81) * t80 + (-t64 * t81 + t65 * t83) * mrSges(5,3)) * t80;
T = t1;
