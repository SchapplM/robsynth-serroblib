% Calculate kinetic energy for
% S5RPRPR3
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
% Datum: 2022-01-23 09:21
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR3_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR3_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2022-01-23 09:20:33
% EndTime: 2022-01-23 09:20:33
% DurationCPUTime: 0.17s
% Computational Cost: add. (186->52), mult. (328->85), div. (0->0), fcn. (160->8), ass. (0->27)
t81 = cos(pkin(8));
t72 = (pkin(1) * t81 + pkin(2)) * qJD(1);
t83 = sin(qJ(3));
t85 = cos(qJ(3));
t79 = sin(pkin(8));
t90 = pkin(1) * qJD(1) * t79;
t70 = t83 * t72 + t85 * t90;
t77 = qJD(1) + qJD(3);
t68 = qJ(4) * t77 + t70;
t78 = sin(pkin(9));
t80 = cos(pkin(9));
t64 = -t80 * qJD(2) + t68 * t78;
t93 = t64 ^ 2;
t92 = m(3) / 0.2e1;
t91 = t77 * t80;
t69 = t72 * t85 - t83 * t90;
t89 = qJD(4) - t69;
t86 = qJD(2) ^ 2;
t84 = cos(qJ(5));
t82 = sin(qJ(5));
t73 = qJD(5) - t91;
t67 = -pkin(3) * t77 + t89;
t66 = qJD(2) * t78 + t68 * t80;
t63 = (-pkin(4) * t80 - pkin(7) * t78 - pkin(3)) * t77 + t89;
t62 = t63 * t82 + t66 * t84;
t61 = t63 * t84 - t66 * t82;
t1 = m(6) * (t61 ^ 2 + t62 ^ 2 + t93) / 0.2e1 + m(4) * (t69 ^ 2 + t70 ^ 2 + t86) / 0.2e1 + t86 * t92 + m(5) * (t66 ^ 2 + t67 ^ 2 + t93) / 0.2e1 + (t61 * mrSges(6,1) - t62 * mrSges(6,2) + Ifges(6,3) * t73 / 0.2e1) * t73 + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t81 * mrSges(3,1) - t79 * mrSges(3,2) + (t79 ^ 2 + t81 ^ 2) * t92 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (-t70 * mrSges(4,2) + t69 * mrSges(4,1) + Ifges(4,3) * t77 / 0.2e1 + (-t67 * mrSges(5,1) + t66 * mrSges(5,3) + Ifges(5,2) * t91 / 0.2e1) * t80 + (t67 * mrSges(5,2) + (Ifges(5,4) * t80 + (Ifges(6,1) * t84 ^ 2 / 0.2e1 + Ifges(5,1) / 0.2e1 + (-Ifges(6,4) * t84 + Ifges(6,2) * t82 / 0.2e1) * t82) * t78) * t77 + (mrSges(6,1) * t82 + mrSges(6,2) * t84 + mrSges(5,3)) * t64 + (-t61 * t84 - t62 * t82) * mrSges(6,3) + t73 * (Ifges(6,5) * t84 - Ifges(6,6) * t82)) * t78) * t77;
T = t1;
