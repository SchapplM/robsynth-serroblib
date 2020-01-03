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
% Datum: 2020-01-03 11:54
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
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
% StartTime: 2020-01-03 11:53:44
% EndTime: 2020-01-03 11:53:45
% DurationCPUTime: 0.20s
% Computational Cost: add. (217->61), mult. (369->100), div. (0->0), fcn. (186->8), ass. (0->30)
t98 = m(3) / 0.2e1;
t84 = qJD(1) + qJD(3);
t91 = cos(qJ(4));
t97 = t84 * t91;
t86 = cos(pkin(9));
t78 = (pkin(1) * t86 + pkin(2)) * qJD(1);
t89 = sin(qJ(3));
t92 = cos(qJ(3));
t85 = sin(pkin(9));
t96 = pkin(1) * qJD(1) * t85;
t74 = t89 * t78 + t92 * t96;
t72 = t84 * pkin(7) + t74;
t88 = sin(qJ(4));
t68 = t88 * qJD(2) + t91 * t72;
t73 = t92 * t78 - t89 * t96;
t93 = qJD(2) ^ 2;
t90 = cos(qJ(5));
t87 = sin(qJ(5));
t83 = qJD(4) + qJD(5);
t82 = t91 * qJD(2);
t76 = (t87 * t91 + t88 * t90) * t84;
t75 = (-t87 * t88 + t90 * t91) * t84;
t71 = -t84 * pkin(3) - t73;
t69 = (-pkin(4) * t91 - pkin(3)) * t84 - t73;
t67 = -t88 * t72 + t82;
t66 = pkin(8) * t97 + t68;
t65 = qJD(4) * pkin(4) + t82 + (-pkin(8) * t84 - t72) * t88;
t64 = t87 * t65 + t90 * t66;
t63 = t90 * t65 - t87 * t66;
t1 = m(5) * (t67 ^ 2 + t68 ^ 2 + t71 ^ 2) / 0.2e1 + m(6) * (t63 ^ 2 + t64 ^ 2 + t69 ^ 2) / 0.2e1 + m(4) * (t73 ^ 2 + t74 ^ 2 + t93) / 0.2e1 + t93 * t98 + (t63 * mrSges(6,1) - t64 * mrSges(6,2) + Ifges(6,3) * t83 / 0.2e1) * t83 + (t67 * mrSges(5,1) - t68 * mrSges(5,2) + Ifges(5,3) * qJD(4) / 0.2e1) * qJD(4) + (t69 * mrSges(6,2) - t63 * mrSges(6,3) + Ifges(6,5) * t83 + Ifges(6,1) * t76 / 0.2e1) * t76 + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t86 * mrSges(3,1) - t85 * mrSges(3,2) + (t85 ^ 2 + t86 ^ 2) * t98 * pkin(1)) * pkin(1)) * qJD(1) ^ 2 + (-t69 * mrSges(6,1) + t64 * mrSges(6,3) + Ifges(6,4) * t76 + Ifges(6,6) * t83 + Ifges(6,2) * t75 / 0.2e1) * t75 + (-t74 * mrSges(4,2) + t73 * mrSges(4,1) + Ifges(4,3) * t84 / 0.2e1 + (-t71 * mrSges(5,1) + t68 * mrSges(5,3) + Ifges(5,6) * qJD(4) + Ifges(5,2) * t97 / 0.2e1) * t91 + (t71 * mrSges(5,2) - t67 * mrSges(5,3) + Ifges(5,5) * qJD(4) + (Ifges(5,4) * t91 + Ifges(5,1) * t88 / 0.2e1) * t84) * t88) * t84;
T = t1;
