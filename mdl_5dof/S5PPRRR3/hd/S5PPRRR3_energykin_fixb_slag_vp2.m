% Calculate kinetic energy for
% S5PPRRR3
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
% Datum: 2019-12-05 15:17
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5PPRRR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5PPRRR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 15:16:18
% EndTime: 2019-12-05 15:16:18
% DurationCPUTime: 0.16s
% Computational Cost: add. (129->56), mult. (288->94), div. (0->0), fcn. (177->8), ass. (0->29)
t86 = sin(qJ(3));
t89 = cos(qJ(3));
t82 = sin(pkin(9));
t94 = qJD(1) * t82;
t76 = t86 * qJD(2) + t89 * t94;
t83 = cos(pkin(9));
t93 = qJD(1) * t83;
t88 = cos(qJ(4));
t92 = t88 * qJD(3);
t91 = t88 * t93;
t75 = qJD(2) * t89 - t86 * t94;
t72 = qJD(3) * pkin(6) + t76;
t85 = sin(qJ(4));
t68 = t88 * t72 - t85 * t93;
t90 = qJD(1) ^ 2;
t87 = cos(qJ(5));
t84 = sin(qJ(5));
t81 = qJD(4) + qJD(5);
t79 = t83 ^ 2 * t90;
t74 = (t84 * t88 + t85 * t87) * qJD(3);
t73 = (-t84 * t85 + t87 * t88) * qJD(3);
t71 = -qJD(3) * pkin(3) - t75;
t69 = (-pkin(4) * t88 - pkin(3)) * qJD(3) - t75;
t67 = -t72 * t85 - t91;
t66 = pkin(7) * t92 + t68;
t65 = -t91 + qJD(4) * pkin(4) + (-pkin(7) * qJD(3) - t72) * t85;
t64 = t65 * t84 + t66 * t87;
t63 = t65 * t87 - t66 * t84;
t1 = m(3) * (t82 ^ 2 * t90 + qJD(2) ^ 2 + t79) / 0.2e1 + m(2) * t90 / 0.2e1 + m(5) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t79) / 0.2e1 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t74 / 0.2e1) * t74 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t74 + Ifges(6,6) * t81 + Ifges(6,2) * t73 / 0.2e1) * t73 + (t75 * mrSges(4,1) - t76 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1 + (-t71 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t92 / 0.2e1) * t88 + (t71 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t88 + Ifges(5,1) * t85 / 0.2e1) * qJD(3)) * t85) * qJD(3);
T = t1;
