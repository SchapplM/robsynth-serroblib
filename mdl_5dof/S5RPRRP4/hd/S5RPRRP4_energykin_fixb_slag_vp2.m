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
% Datum: 2019-12-05 18:07
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
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
% StartTime: 2019-12-05 18:05:25
% EndTime: 2019-12-05 18:05:25
% DurationCPUTime: 0.40s
% Computational Cost: add. (265->77), mult. (652->121), div. (0->0), fcn. (406->6), ass. (0->34)
t89 = qJD(1) ^ 2;
t96 = t89 * qJ(2) ^ 2;
t83 = sin(pkin(8));
t84 = cos(pkin(8));
t72 = qJD(2) + (-pkin(2) * t84 - pkin(6) * t83 - pkin(1)) * qJD(1);
t88 = cos(qJ(3));
t71 = t88 * t72;
t94 = qJD(1) * t84;
t77 = qJD(3) - t94;
t86 = sin(qJ(3));
t62 = pkin(3) * t77 + t71 + (-pkin(7) * t83 * t88 - qJ(2) * t84 * t86) * qJD(1);
t93 = qJD(1) * qJ(2);
t91 = t84 * t93;
t67 = t86 * t72 + t88 * t91;
t95 = qJD(1) * t83;
t92 = t86 * t95;
t65 = -pkin(7) * t92 + t67;
t85 = sin(qJ(4));
t87 = cos(qJ(4));
t59 = t85 * t62 + t87 * t65;
t73 = pkin(3) * t92 + t83 * t93;
t58 = t87 * t62 - t65 * t85;
t82 = t84 ^ 2;
t81 = t83 ^ 2;
t80 = -qJD(1) * pkin(1) + qJD(2);
t78 = t81 * t96;
t75 = qJD(4) + t77;
t69 = (-t85 * t86 + t87 * t88) * t95;
t68 = (-t85 * t88 - t86 * t87) * t95;
t66 = -t86 * t91 + t71;
t64 = -pkin(4) * t68 + qJD(5) + t73;
t57 = qJ(5) * t68 + t59;
t56 = pkin(4) * t75 - qJ(5) * t69 + t58;
t1 = m(6) * (t56 ^ 2 + t57 ^ 2 + t64 ^ 2) / 0.2e1 + m(4) * (t66 ^ 2 + t67 ^ 2 + t78) / 0.2e1 + m(3) * (t80 ^ 2 + t82 * t96 + t78) / 0.2e1 + m(5) * (t58 ^ 2 + t59 ^ 2 + t73 ^ 2) / 0.2e1 + (t66 * mrSges(4,1) - t67 * mrSges(4,2) + Ifges(4,3) * t77 / 0.2e1) * t77 + (t58 * mrSges(5,1) + t56 * mrSges(6,1) - t59 * mrSges(5,2) - t57 * mrSges(6,2) + (Ifges(5,3) / 0.2e1 + Ifges(6,3) / 0.2e1) * t75) * t75 + (t73 * mrSges(5,2) + t64 * mrSges(6,2) - t58 * mrSges(5,3) - t56 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t69 + (Ifges(5,5) + Ifges(6,5)) * t75) * t69 + (-t73 * mrSges(5,1) - t64 * mrSges(6,1) + t59 * mrSges(5,3) + t57 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t68 + (Ifges(5,6) + Ifges(6,6)) * t75 + (Ifges(5,4) + Ifges(6,4)) * t69) * t68 + ((-t80 * mrSges(3,1) + Ifges(3,2) * t94 / 0.2e1) * t84 + (t80 * mrSges(3,2) + (Ifges(3,4) * t84 + (Ifges(3,1) / 0.2e1 + (qJ(2) * mrSges(4,2) + Ifges(4,1) * t88 / 0.2e1) * t88 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t88 + Ifges(4,2) * t86 / 0.2e1) * t86) * t83) * qJD(1) + (-t66 * t88 - t67 * t86) * mrSges(4,3) + t77 * (Ifges(4,5) * t88 - Ifges(4,6) * t86)) * t83) * qJD(1) + (Ifges(2,3) / 0.2e1 + (t81 + t82) * qJ(2) * mrSges(3,3)) * t89;
T = t1;
