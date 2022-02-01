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
% Datum: 2022-01-20 09:52
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-20 09:51:16
% EndTime: 2022-01-20 09:51:17
% DurationCPUTime: 0.20s
% Computational Cost: add. (224->54), mult. (353->93), div. (0->0), fcn. (184->8), ass. (0->29)
t83 = qJD(1) + qJD(2);
t97 = pkin(7) * t83;
t91 = cos(qJ(2));
t96 = pkin(1) * qJD(1);
t77 = t83 * pkin(2) + t91 * t96;
t85 = sin(pkin(8));
t87 = cos(pkin(8));
t89 = sin(qJ(2));
t95 = t89 * t96;
t73 = t85 * t77 + t87 * t95;
t71 = t83 * qJ(4) + t73;
t84 = sin(pkin(9));
t86 = cos(pkin(9));
t67 = t84 * qJD(3) + t86 * t71;
t72 = t87 * t77 - t85 * t95;
t94 = qJD(4) - t72;
t90 = cos(qJ(5));
t88 = sin(qJ(5));
t81 = t86 * qJD(3);
t75 = (t84 * t90 + t86 * t88) * t83;
t74 = (-t84 * t88 + t86 * t90) * t83;
t70 = -t83 * pkin(3) + t94;
t68 = (-pkin(4) * t86 - pkin(3)) * t83 + t94;
t66 = -t84 * t71 + t81;
t65 = t86 * t97 + t67;
t64 = t81 + (-t71 - t97) * t84;
t63 = t88 * t64 + t90 * t65;
t62 = t90 * t64 - t88 * t65;
t1 = m(4) * (qJD(3) ^ 2 + t72 ^ 2 + t73 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t70 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + (Ifges(2,3) / 0.2e1 + m(3) * (t89 ^ 2 + t91 ^ 2) * pkin(1) ^ 2 / 0.2e1) * qJD(1) ^ 2 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,1) * t75 / 0.2e1) * t75 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t75 + Ifges(6,2) * t74 / 0.2e1) * t74 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,5) * t75 + Ifges(6,6) * t74 + Ifges(6,3) * qJD(5) / 0.2e1) * qJD(5) + (-t73 * mrSges(4,2) + t72 * mrSges(4,1) + t70 * (-mrSges(5,1) * t86 + mrSges(5,2) * t84) + (mrSges(3,1) * t91 - mrSges(3,2) * t89) * t96 + (-t66 * t84 + t67 * t86) * mrSges(5,3) + (Ifges(3,3) / 0.2e1 + Ifges(4,3) / 0.2e1 + Ifges(5,2) * t86 ^ 2 / 0.2e1 + (Ifges(5,4) * t86 + Ifges(5,1) * t84 / 0.2e1) * t84) * t83) * t83;
T = t1;
