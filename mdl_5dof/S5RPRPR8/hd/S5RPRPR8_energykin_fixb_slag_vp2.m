% Calculate kinetic energy for
% S5RPRPR8
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
% Datum: 2019-12-31 18:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR8_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR8_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:21:14
% EndTime: 2019-12-31 18:21:15
% DurationCPUTime: 0.35s
% Computational Cost: add. (252->73), mult. (556->118), div. (0->0), fcn. (324->8), ass. (0->32)
t100 = m(3) / 0.2e1;
t88 = sin(pkin(8));
t83 = (pkin(1) * t88 + pkin(6)) * qJD(1);
t92 = sin(qJ(3));
t94 = cos(qJ(3));
t78 = t92 * qJD(2) + t94 * t83;
t75 = qJD(3) * qJ(4) + t78;
t90 = cos(pkin(8));
t97 = -pkin(1) * t90 - pkin(2);
t76 = (-pkin(3) * t94 - qJ(4) * t92 + t97) * qJD(1);
t87 = sin(pkin(9));
t89 = cos(pkin(9));
t67 = t89 * t75 + t87 * t76;
t99 = qJD(1) * t92;
t98 = qJD(1) * t94;
t66 = -t75 * t87 + t89 * t76;
t77 = qJD(2) * t94 - t92 * t83;
t74 = -qJD(3) * pkin(3) + qJD(4) - t77;
t93 = cos(qJ(5));
t91 = sin(qJ(5));
t85 = qJD(5) - t98;
t84 = t97 * qJD(1);
t80 = qJD(3) * t87 + t89 * t99;
t79 = qJD(3) * t89 - t87 * t99;
t70 = t79 * t91 + t80 * t93;
t69 = t79 * t93 - t80 * t91;
t68 = -pkin(4) * t79 + t74;
t65 = pkin(7) * t79 + t67;
t64 = -pkin(4) * t98 - pkin(7) * t80 + t66;
t63 = t64 * t91 + t65 * t93;
t62 = t64 * t93 - t65 * t91;
t1 = qJD(2) ^ 2 * t100 + m(4) * (t77 ^ 2 + t78 ^ 2 + t84 ^ 2) / 0.2e1 + m(5) * (t66 ^ 2 + t67 ^ 2 + t74 ^ 2) / 0.2e1 + m(6) * (t62 ^ 2 + t63 ^ 2 + t68 ^ 2) / 0.2e1 + (t62 * mrSges(6,1) - t63 * mrSges(6,2) + Ifges(6,3) * t85 / 0.2e1) * t85 + (t74 * mrSges(5,2) - t66 * mrSges(5,3) + Ifges(5,1) * t80 / 0.2e1) * t80 + (-t74 * mrSges(5,1) + t67 * mrSges(5,3) + Ifges(5,4) * t80 + Ifges(5,2) * t79 / 0.2e1) * t79 + (t68 * mrSges(6,2) - t62 * mrSges(6,3) + Ifges(6,5) * t85 + Ifges(6,1) * t70 / 0.2e1) * t70 + (-t68 * mrSges(6,1) + t63 * mrSges(6,3) + Ifges(6,4) * t70 + Ifges(6,6) * t85 + Ifges(6,2) * t69 / 0.2e1) * t69 + (t77 * mrSges(4,1) - t78 * mrSges(4,2) + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + ((t84 * mrSges(4,2) - t77 * mrSges(4,3) + Ifges(4,1) * t99 / 0.2e1) * t92 + (-t84 * mrSges(4,1) - t66 * mrSges(5,1) + t67 * mrSges(5,2) + t78 * mrSges(4,3) - Ifges(5,5) * t80 - Ifges(5,6) * t79) * t94 + qJD(3) * (Ifges(4,5) * t92 + Ifges(4,6) * t94) + (Ifges(2,3) / 0.2e1 + Ifges(3,3) / 0.2e1 + (t90 * mrSges(3,1) - t88 * mrSges(3,2) + (t88 ^ 2 + t90 ^ 2) * t100 * pkin(1)) * pkin(1) + (Ifges(4,4) * t92 + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t94) * t94) * qJD(1)) * qJD(1);
T = t1;
