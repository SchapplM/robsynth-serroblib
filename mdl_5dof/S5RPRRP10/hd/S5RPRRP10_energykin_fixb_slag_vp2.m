% Calculate kinetic energy for
% S5RPRRP10
% Use Code from Maple symbolic Code Generation
% 
% Input:
% qJ [5x1]
%   Generalized joint coordinates (joint angles)
% qJD [5x1]
%   Generalized joint velocities
% pkin [8x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,d1,d3,d4,theta2]';
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
% Datum: 2019-12-31 18:52
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP10_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP10_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2019-12-31 18:50:51
% EndTime: 2019-12-31 18:50:51
% DurationCPUTime: 0.38s
% Computational Cost: add. (297->76), mult. (724->112), div. (0->0), fcn. (482->6), ass. (0->29)
t84 = cos(pkin(8));
t96 = t84 ^ 2;
t95 = qJD(1) * (pkin(6) + qJ(2));
t83 = sin(pkin(8));
t86 = sin(qJ(3));
t88 = cos(qJ(3));
t75 = (t83 * t86 - t84 * t88) * qJD(1);
t94 = m(3) / 0.2e1;
t76 = (t83 * t88 + t84 * t86) * qJD(1);
t79 = qJD(2) + (-pkin(2) * t84 - pkin(1)) * qJD(1);
t63 = pkin(3) * t75 - pkin(7) * t76 + t79;
t77 = t83 * t95;
t78 = t84 * t95;
t68 = -t86 * t77 + t88 * t78;
t66 = qJD(3) * pkin(7) + t68;
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t59 = t85 * t63 + t87 * t66;
t58 = t87 * t63 - t66 * t85;
t67 = -t77 * t88 - t86 * t78;
t65 = -qJD(3) * pkin(3) - t67;
t80 = -qJD(1) * pkin(1) + qJD(2);
t71 = qJD(4) + t75;
t70 = qJD(3) * t85 + t76 * t87;
t69 = qJD(3) * t87 - t76 * t85;
t60 = -pkin(4) * t69 + qJD(5) + t65;
t57 = qJ(5) * t69 + t59;
t56 = pkin(4) * t71 - qJ(5) * t70 + t58;
t1 = t80 ^ 2 * t94 + m(4) * (t67 ^ 2 + t68 ^ 2 + t79 ^ 2) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t65 ^ 2) / 0.2e1 + m(6) * (t56 ^ 2 + t57 ^ 2 + t60 ^ 2) / 0.2e1 + (t79 * mrSges(4,2) - t67 * mrSges(4,3) + Ifges(4,1) * t76 / 0.2e1) * t76 - (-t79 * mrSges(4,1) + t68 * mrSges(4,3) + Ifges(4,4) * t76 - Ifges(4,2) * t75 / 0.2e1) * t75 + (t67 * mrSges(4,1) - t68 * mrSges(4,2) + Ifges(4,5) * t76 - Ifges(4,6) * t75 + Ifges(4,3) * qJD(3) / 0.2e1) * qJD(3) + (t58 * mrSges(5,1) + t56 * mrSges(6,1) - t59 * mrSges(5,2) - t57 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t71) * t71 + (t65 * mrSges(5,2) + t60 * mrSges(6,2) - t58 * mrSges(5,3) - t56 * mrSges(6,3) + (Ifges(6,1) / 0.2e1 + Ifges(5,1) / 0.2e1) * t70 + (Ifges(5,5) + Ifges(6,5)) * t71) * t70 + (-t65 * mrSges(5,1) - t60 * mrSges(6,1) + t59 * mrSges(5,3) + t57 * mrSges(6,3) + (Ifges(6,2) / 0.2e1 + Ifges(5,2) / 0.2e1) * t69 + (Ifges(5,6) + Ifges(6,6)) * t71 + (Ifges(5,4) + Ifges(6,4)) * t70) * t69 + (t80 * (-mrSges(3,1) * t84 + mrSges(3,2) * t83) + (Ifges(2,3) / 0.2e1 + (qJ(2) * t94 + mrSges(3,3)) * (t83 ^ 2 + t96) * qJ(2) + Ifges(3,2) * t96 / 0.2e1 + (Ifges(3,4) * t84 + Ifges(3,1) * t83 / 0.2e1) * t83) * qJD(1)) * qJD(1);
T = t1;
