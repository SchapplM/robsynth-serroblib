% Calculate kinetic energy for
% S5PPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d3,d4,d5,theta1,theta2]';
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
% Datum: 2019-12-05 15:15
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR2_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR2_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:14:07
% EndTime: 2019-12-05 15:14:07
% DurationCPUTime: 0.16s
% Computational Cost: add. (129->56), mult. (294->94), div. (0->0), fcn. (184->8), ass. (0->27)
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t94 = qJD(1) * cos(qJ(3));
t95 = qJD(1) * sin(qJ(3));
t74 = t83 * t94 + t84 * t95;
t72 = qJD(3) * pkin(6) + t74;
t86 = sin(qJ(4));
t89 = cos(qJ(4));
t68 = t86 * qJD(2) + t89 * t72;
t93 = qJD(3) * t89;
t73 = -t83 * t95 + t84 * t94;
t92 = qJD(1) ^ 2;
t91 = qJD(2) ^ 2;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t82 = qJD(4) + qJD(5);
t81 = t89 * qJD(2);
t76 = (t85 * t89 + t86 * t88) * qJD(3);
t75 = (-t85 * t86 + t88 * t89) * qJD(3);
t71 = -qJD(3) * pkin(3) - t73;
t69 = (-pkin(4) * t89 - pkin(3)) * qJD(3) - t73;
t67 = -t72 * t86 + t81;
t66 = pkin(7) * t93 + t68;
t65 = qJD(4) * pkin(4) + t81 + (-pkin(7) * qJD(3) - t72) * t86;
t64 = t65 * t85 + t66 * t88;
t63 = t65 * t88 - t66 * t85;
t1 = m(3) * (t91 + (t83 ^ 2 + t84 ^ 2) * t92) / 0.2e1 + m(2) * t92 / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t91) / 0.2e1 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t82 / 0.2e1) * t82 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t82 + Ifges(6,1) * t76 / 0.2e1) * t76 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t76 + Ifges(6,6) * t82 + Ifges(6,2) * t75 / 0.2e1) * t75 + (t73 * mrSges(4,1) - t74 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t71 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t93 / 0.2e1) * t89 + (t71 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t89 + Ifges(5,1) * t86 / 0.2e1) * qJD(3)) * t86) * qJD(3);
T = t1;
