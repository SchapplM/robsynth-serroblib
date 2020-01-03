% Calculate kinetic energy for
% S5RPRPR13
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d5,theta2]';
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
% Datum: 2019-12-31 18:33
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRPR13_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRPR13_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:31:52
% EndTime: 2019-12-31 18:31:53
% DurationCPUTime: 0.37s
% Computational Cost: add. (253->77), mult. (622->112), div. (0->0), fcn. (392->6), ass. (0->34)
t83 = cos(pkin(8));
t98 = t83 ^ 2;
t97 = m(3) / 0.2e1;
t96 = pkin(3) + pkin(7);
t95 = cos(qJ(3));
t94 = pkin(6) + qJ(2);
t82 = sin(pkin(8));
t92 = qJD(1) * t82;
t75 = t94 * t92;
t91 = qJD(1) * t83;
t76 = t94 * t91;
t85 = sin(qJ(3));
t66 = -t85 * t75 + t95 * t76;
t65 = -t75 * t95 - t85 * t76;
t64 = -qJD(3) * qJ(4) - t66;
t90 = qJD(4) - t65;
t77 = qJD(2) + (-pkin(2) * t83 - pkin(1)) * qJD(1);
t74 = (t82 * t95 + t83 * t85) * qJD(1);
t89 = -qJ(4) * t74 + t77;
t86 = cos(qJ(5));
t84 = sin(qJ(5));
t79 = -qJD(1) * pkin(1) + qJD(2);
t73 = t85 * t92 - t91 * t95;
t69 = qJD(5) + t74;
t68 = qJD(3) * t86 + t73 * t84;
t67 = -qJD(3) * t84 + t73 * t86;
t63 = -qJD(3) * pkin(3) + t90;
t62 = pkin(3) * t73 + t89;
t61 = -pkin(4) * t73 - t64;
t60 = t74 * pkin(4) - qJD(3) * t96 + t90;
t59 = t73 * t96 + t89;
t58 = t59 * t86 + t60 * t84;
t57 = -t59 * t84 + t60 * t86;
t1 = t79 ^ 2 * t97 + m(4) * (t65 ^ 2 + t66 ^ 2 + t77 ^ 2) / 0.2e1 + m(6) * (t57 ^ 2 + t58 ^ 2 + t61 ^ 2) / 0.2e1 + m(5) * (t62 ^ 2 + t63 ^ 2 + t64 ^ 2) / 0.2e1 + (t57 * mrSges(6,1) - t58 * mrSges(6,2) + Ifges(6,3) * t69 / 0.2e1) * t69 + (t61 * mrSges(6,2) - t57 * mrSges(6,3) + Ifges(6,5) * t69 + Ifges(6,1) * t68 / 0.2e1) * t68 + (-t61 * mrSges(6,1) + t58 * mrSges(6,3) + Ifges(6,4) * t68 + Ifges(6,6) * t69 + Ifges(6,2) * t67 / 0.2e1) * t67 + (t63 * mrSges(5,1) + t77 * mrSges(4,2) - t65 * mrSges(4,3) - t62 * mrSges(5,3) + (Ifges(4,1) / 0.2e1 + Ifges(5,2) / 0.2e1) * t74) * t74 + (t77 * mrSges(4,1) + t64 * mrSges(5,1) - t62 * mrSges(5,2) - t66 * mrSges(4,3) + (Ifges(4,2) / 0.2e1 + Ifges(5,3) / 0.2e1) * t73 + (-Ifges(4,4) - Ifges(5,6)) * t74) * t73 + (t65 * mrSges(4,1) - t66 * mrSges(4,2) + t63 * mrSges(5,2) - t64 * mrSges(5,3) + (Ifges(4,3) / 0.2e1 + Ifges(5,1) / 0.2e1) * qJD(3) + (-Ifges(5,4) + Ifges(4,5)) * t74 + (Ifges(5,5) - Ifges(4,6)) * t73) * qJD(3) + (t79 * (-mrSges(3,1) * t83 + mrSges(3,2) * t82) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t97 + mrSges(3,3)) * (t82 ^ 2 + t98) * qJ(2) + Ifges(3,2) * t98 / 0.2e1 + (Ifges(3,4) * t83 + Ifges(3,1) * t82 / 0.2e1) * t82) * qJD(1)) * qJD(1);
T = t1;
