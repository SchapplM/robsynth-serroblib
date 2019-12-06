% Calculate kinetic energy for
% S5RRRPR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [9x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d2,d3,d5,theta4]';
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
% Datum: 2019-12-05 18:43
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RRRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RRRPR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-05 18:42:15
% EndTime: 2019-12-05 18:42:16
% DurationCPUTime: 0.25s
% Computational Cost: add. (400->71), mult. (579->119), div. (0->0), fcn. (346->8), ass. (0->32)
t82 = qJD(1) + qJD(2);
t97 = t82 / 0.2e1;
t87 = sin(qJ(2));
t95 = qJD(1) * pkin(1);
t79 = t82 * pkin(7) + t87 * t95;
t96 = t79 * mrSges(4,3);
t86 = sin(qJ(3));
t93 = qJ(4) * t82 + t79;
t73 = qJD(3) * pkin(3) - t93 * t86;
t89 = cos(qJ(3));
t74 = t93 * t89;
t83 = sin(pkin(9));
t84 = cos(pkin(9));
t66 = t83 * t73 + t84 * t74;
t90 = cos(qJ(2));
t94 = t90 * t95;
t65 = t84 * t73 - t83 * t74;
t75 = -t94 + qJD(4) + (-pkin(3) * t89 - pkin(2)) * t82;
t88 = cos(qJ(5));
t85 = sin(qJ(5));
t81 = qJD(3) + qJD(5);
t80 = -t82 * pkin(2) - t94;
t77 = (t83 * t89 + t84 * t86) * t82;
t76 = (-t83 * t86 + t84 * t89) * t82;
t69 = -t76 * pkin(4) + t75;
t68 = t85 * t76 + t88 * t77;
t67 = t88 * t76 - t85 * t77;
t64 = t76 * pkin(8) + t66;
t63 = qJD(3) * pkin(4) - t77 * pkin(8) + t65;
t62 = t85 * t63 + t88 * t64;
t61 = t88 * t63 - t85 * t64;
t1 = m(4) * (t80 ^ 2 + (t86 ^ 2 + t89 ^ 2) * t79 ^ 2) / 0.2e1 + m(6) * (t61 ^ 2 + t62 ^ 2 + t69 ^ 2) / 0.2e1 + m(5) * (t65 ^ 2 + t66 ^ 2 + t75 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t87 ^ 2 + t90 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t81 / 0.2e1) * t81 + (t75 * mrSges(5,2) - t65 * mrSges(5,3) + Ifges(5,1) * t77 / 0.2e1) * t77 + (-t75 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,4) * t77 + Ifges(5,2) * t76 / 0.2e1) * t76 + (t69 * mrSges(6,2) - t61 * mrSges(6,3) + Ifges(6,5) * t81 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t69 * mrSges(6,1) + t62 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t81 + Ifges(6,2) * t67 / 0.2e1) * t67 + (Ifges(3,3) * t97 + (mrSges(3,1) * t90 - mrSges(3,2) * t87) * t95 + (-t80 * mrSges(4,1) + (Ifges(4,2) * t97 + t96) * t89) * t89 + (Ifges(4,4) * t89 * t82 + t80 * mrSges(4,2) + (Ifges(4,1) * t97 + t96) * t86) * t86) * t82 + (t65 * mrSges(5,1) - t66 * mrSges(5,2) + Ifges(5,5) * t77 + Ifges(5,6) * t76 + (Ifges(4,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * qJD(3) + (Ifges(4,5) * t86 + Ifges(4,6) * t89) * t82 + (-mrSges(4,1) * t86 - mrSges(4,2) * t89) * t79) * qJD(3);
T = t1;
