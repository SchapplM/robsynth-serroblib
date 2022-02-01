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
% m [6x1]
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
% Datum: 2022-01-20 09:13
% Revision: 008671b0a00594318b890887636eaaff83cd5e2f (2021-12-12)
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
% StartTime: 2022-01-20 09:12:09
% EndTime: 2022-01-20 09:12:09
% DurationCPUTime: 0.29s
% Computational Cost: add. (197->66), mult. (483->110), div. (0->0), fcn. (276->8), ass. (0->30)
t96 = m(3) / 0.2e1;
t83 = sin(pkin(8));
t86 = cos(pkin(8));
t87 = cos(pkin(7));
t92 = -pkin(1) * t87 - pkin(2);
t71 = qJD(3) + (-pkin(3) * t86 - qJ(4) * t83 + t92) * qJD(1);
t84 = sin(pkin(7));
t80 = (pkin(1) * t84 + qJ(3)) * qJD(1);
t77 = qJD(2) * t83 + t80 * t86;
t82 = sin(pkin(9));
t85 = cos(pkin(9));
t67 = t82 * t71 + t85 * t77;
t95 = qJD(1) * t83;
t94 = qJD(1) * t86;
t93 = t82 * t95;
t66 = t85 * t71 - t77 * t82;
t76 = qJD(2) * t86 - t83 * t80;
t75 = qJD(4) - t76;
t89 = cos(qJ(5));
t88 = sin(qJ(5));
t81 = qJD(5) - t94;
t79 = t92 * qJD(1) + qJD(3);
t74 = (-t82 * t88 + t85 * t89) * t95;
t73 = (-t82 * t89 - t85 * t88) * t95;
t68 = pkin(4) * t93 + t75;
t65 = -pkin(6) * t93 + t67;
t64 = (-pkin(6) * t83 * t85 - pkin(4) * t86) * qJD(1) + t66;
t63 = t64 * t88 + t65 * t89;
t62 = t64 * t89 - t65 * t88;
t1 = qJD(2) ^ 2 * t96 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t75 ^ 2) / 0.2e1 + m(4) * (t76 ^ 2 + t77 ^ 2 + t79 ^ 2) / 0.2e1 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,6) * t81 + Ifges(6,2) * t73 / 0.2e1) * t73 + ((-t76 * mrSges(4,3) + t75 * (mrSges(5,1) * t82 + mrSges(5,2) * t85) + t79 * mrSges(4,2) + (Ifges(5,1) * t85 ^ 2 / 0.2e1 + Ifges(4,1) / 0.2e1 + (-Ifges(5,4) * t85 + Ifges(5,2) * t82 / 0.2e1) * t82) * t95 + (-t66 * t85 - t67 * t82) * mrSges(5,3)) * t83 + (t77 * mrSges(4,3) + t67 * mrSges(5,2) - t79 * mrSges(4,1) - t66 * mrSges(5,1) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t94 + (-Ifges(5,5) * t85 + Ifges(5,6) * t82 + Ifges(4,4)) * t95) * t86 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t87 * mrSges(3,1) - t84 * mrSges(3,2) + (t84 ^ 2 + t87 ^ 2) * t96 * pkin(1)) * pkin(1)) * qJD(1)) * qJD(1);
T = t1;
