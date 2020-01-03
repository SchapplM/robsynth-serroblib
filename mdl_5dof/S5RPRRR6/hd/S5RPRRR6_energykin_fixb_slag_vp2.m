% Calculate kinetic energy for
% S5RPRRR6
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
% Datum: 2019-12-31 19:02
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRR6_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR6_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 19:00:54
% EndTime: 2019-12-31 19:00:54
% DurationCPUTime: 0.38s
% Computational Cost: add. (260->74), mult. (550->122), div. (0->0), fcn. (326->8), ass. (0->33)
t90 = sin(qJ(4));
t91 = sin(qJ(3));
t93 = cos(qJ(4));
t94 = cos(qJ(3));
t78 = (t90 * t91 - t93 * t94) * qJD(1);
t100 = m(3) / 0.2e1;
t87 = sin(pkin(9));
t82 = (pkin(1) * t87 + pkin(6)) * qJD(1);
t85 = t94 * qJD(2);
t99 = pkin(7) * qJD(1);
t71 = qJD(3) * pkin(3) + t85 + (-t82 - t99) * t91;
t76 = t91 * qJD(2) + t94 * t82;
t74 = t94 * t99 + t76;
t67 = t90 * t71 + t93 * t74;
t88 = cos(pkin(9));
t98 = -pkin(1) * t88 - pkin(2);
t66 = t71 * t93 - t74 * t90;
t80 = (-pkin(3) * t94 + t98) * qJD(1);
t92 = cos(qJ(5));
t89 = sin(qJ(5));
t86 = qJD(3) + qJD(4);
t83 = t98 * qJD(1);
t79 = (t90 * t94 + t91 * t93) * qJD(1);
t77 = qJD(5) + t78;
t75 = -t82 * t91 + t85;
t73 = t79 * t92 + t86 * t89;
t72 = -t79 * t89 + t86 * t92;
t68 = pkin(4) * t78 - pkin(8) * t79 + t80;
t65 = pkin(8) * t86 + t67;
t64 = -pkin(4) * t86 - t66;
t63 = t65 * t92 + t68 * t89;
t62 = -t65 * t89 + t68 * t92;
t1 = qJD(2) ^ 2 * t100 + m(5) * (t66 ^ 2 + t67 ^ 2 + t80 ^ 2) / 0.2e1 + m(4) * (t75 ^ 2 + t76 ^ 2 + t83 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t66 * mrSges(5,1) - t67 * mrSges(5,2) + Ifges(5,3) * t86 / 0.2e1) * t86 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t77 / 0.2e1) * t77 + (t80 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,5) * t86 + Ifges(5,1) * t79 / 0.2e1) * t79 + (t64 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t77 + Ifges(6,1) * t73 / 0.2e1) * t73 - (-t80 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t79 + Ifges(5,6) * t86 - Ifges(5,2) * t78 / 0.2e1) * t78 + (-t64 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t73 + Ifges(6,6) * t77 + Ifges(6,2) * t72 / 0.2e1) * t72 + (t75 * mrSges(4,1) - t76 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t83 * (-mrSges(4,1) * t94 + mrSges(4,2) * t91) + (-t75 * t91 + t76 * t94) * mrSges(4,3) + qJD(3) * (Ifges(4,5) * t91 + Ifges(4,6) * t94) + (Ifges(3,3) / 0.2e1 + Ifges(2,3) / 0.2e1 + (t88 * mrSges(3,1) - t87 * mrSges(3,2) + (t87 ^ 2 + t88 ^ 2) * t100 * pkin(1)) * pkin(1) + Ifges(4,2) * t94 ^ 2 / 0.2e1 + (Ifges(4,4) * t94 + Ifges(4,1) * t91 / 0.2e1) * t91) * qJD(1)) * qJD(1);
T = t1;
