% Calculate kinetic energy for
% S5RPRRP4
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
% Datum: 2020-01-03 11:51
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function T = S5RPRRP4_energykin_fixb_slag_vp2(qJ, qJD, ...
  pkin, m, mrSges, Ifges)
%% Coder Information
%#codegen
%$cgargs {zeros(5,1),zeros(5,1),zeros(8,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(qJ) && all(size(qJ) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: qJ has to be [5x1] (double)');
assert(isreal(qJD) && all(size(qJD) == [5 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: qJD has to be [5x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [8 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: pkin has to be [8x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRP4_energykin_fixb_slag_vp2: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From energy_kinetic_fixb_linkframe_par2_matlab.m
% OptimizationMode: 2
% StartTime: 2020-01-03 11:49:09
% EndTime: 2020-01-03 11:49:09
% DurationCPUTime: 0.46s
% Computational Cost: add. (265->77), mult. (652->121), div. (0->0), fcn. (406->6), ass. (0->34)
t87 = qJD(1) ^ 2;
t94 = t87 * qJ(2) ^ 2;
t81 = sin(pkin(8));
t82 = cos(pkin(8));
t70 = qJD(2) + (-pkin(2) * t82 - pkin(6) * t81 - pkin(1)) * qJD(1);
t86 = cos(qJ(3));
t69 = t86 * t70;
t92 = qJD(1) * t82;
t75 = qJD(3) - t92;
t84 = sin(qJ(3));
t60 = pkin(3) * t75 + t69 + (-pkin(7) * t81 * t86 - qJ(2) * t82 * t84) * qJD(1);
t91 = qJD(1) * qJ(2);
t89 = t82 * t91;
t65 = t84 * t70 + t86 * t89;
t93 = qJD(1) * t81;
t90 = t84 * t93;
t63 = -pkin(7) * t90 + t65;
t83 = sin(qJ(4));
t85 = cos(qJ(4));
t57 = t83 * t60 + t85 * t63;
t71 = pkin(3) * t90 + t81 * t91;
t56 = t85 * t60 - t63 * t83;
t80 = t82 ^ 2;
t79 = t81 ^ 2;
t78 = -qJD(1) * pkin(1) + qJD(2);
t76 = t79 * t94;
t73 = qJD(4) + t75;
t67 = (-t83 * t84 + t85 * t86) * t93;
t66 = (-t83 * t86 - t84 * t85) * t93;
t64 = -t84 * t89 + t69;
t62 = -pkin(4) * t66 + qJD(5) + t71;
t55 = qJ(5) * t66 + t57;
t54 = pkin(4) * t73 - qJ(5) * t67 + t56;
t1 = m(3) * (t78 ^ 2 + t80 * t94 + t76) / 0.2e1 + m(6) * (t54 ^ 2 + t55 ^ 2 + t62 ^ 2) / 0.2e1 + m(4) * (t64 ^ 2 + t65 ^ 2 + t76) / 0.2e1 + m(5) * (t56 ^ 2 + t57 ^ 2 + t71 ^ 2) / 0.2e1 + (t64 * mrSges(4,1) - t65 * mrSges(4,2) + Ifges(4,3) * t75 / 0.2e1) * t75 + (t56 * mrSges(5,1) + t54 * mrSges(6,1) - t57 * mrSges(5,2) - t55 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t73) * t73 + (t71 * mrSges(5,2) + t62 * mrSges(6,2) - t56 * mrSges(5,3) - t54 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t67 + (Ifges(5,5) + Ifges(6,5)) * t73) * t67 + (-t71 * mrSges(5,1) - t62 * mrSges(6,1) + t57 * mrSges(5,3) + t55 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t66 + (Ifges(5,6) + Ifges(6,6)) * t73 + (Ifges(5,4) + Ifges(6,4)) * t67) * t66 + ((-t78 * mrSges(3,1) + Ifges(3,2) * t92 / 0.2e1) * t82 + (t78 * mrSges(3,2) + (Ifges(3,4) * t82 + (Ifges(3,1) / 0.2e1 + (qJ(2) * mrSges(4,2) + Ifges(4,1) * t86 / 0.2e1) * t86 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t86 + Ifges(4,2) * t84 / 0.2e1) * t84) * t81) * qJD(1) + (-t64 * t86 - t65 * t84) * mrSges(4,3) + t75 * (Ifges(4,5) * t86 - Ifges(4,6) * t84)) * t81) * qJD(1) + (Ifges(2,3) / 0.2e1 + (t79 + t80) * qJ(2) * mrSges(3,3)) * t87;
T = t1;
