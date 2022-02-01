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
% Datum: 2022-01-23 09:33
% Revision: fd3771346c4aea32fdeb66112c511235427c26a7 (2022-01-20)
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
% StartTime: 2022-01-23 09:31:43
% EndTime: 2022-01-23 09:31:44
% DurationCPUTime: 0.40s
% Computational Cost: add. (265->77), mult. (652->121), div. (0->0), fcn. (406->6), ass. (0->34)
t90 = qJD(1) ^ 2;
t97 = t90 * qJ(2) ^ 2;
t84 = sin(pkin(8));
t85 = cos(pkin(8));
t73 = qJD(2) + (-pkin(2) * t85 - pkin(6) * t84 - pkin(1)) * qJD(1);
t89 = cos(qJ(3));
t72 = t89 * t73;
t95 = qJD(1) * t85;
t78 = qJD(3) - t95;
t87 = sin(qJ(3));
t63 = pkin(3) * t78 + t72 + (-pkin(7) * t84 * t89 - qJ(2) * t85 * t87) * qJD(1);
t94 = qJ(2) * qJD(1);
t92 = t85 * t94;
t68 = t87 * t73 + t89 * t92;
t96 = qJD(1) * t84;
t93 = t87 * t96;
t66 = -pkin(7) * t93 + t68;
t86 = sin(qJ(4));
t88 = cos(qJ(4));
t60 = t86 * t63 + t88 * t66;
t74 = pkin(3) * t93 + t84 * t94;
t59 = t88 * t63 - t66 * t86;
t83 = t85 ^ 2;
t82 = t84 ^ 2;
t81 = -qJD(1) * pkin(1) + qJD(2);
t79 = t82 * t97;
t76 = qJD(4) + t78;
t70 = (-t86 * t87 + t88 * t89) * t96;
t69 = (-t86 * t89 - t87 * t88) * t96;
t67 = -t87 * t92 + t72;
t65 = -pkin(4) * t69 + qJD(5) + t74;
t58 = qJ(5) * t69 + t60;
t57 = pkin(4) * t76 - qJ(5) * t70 + t59;
t1 = m(6) * (t57 ^ 2 + t58 ^ 2 + t65 ^ 2) / 0.2e1 + m(5) * (t59 ^ 2 + t60 ^ 2 + t74 ^ 2) / 0.2e1 + m(4) * (t67 ^ 2 + t68 ^ 2 + t79) / 0.2e1 + m(3) * (t81 ^ 2 + t83 * t97 + t79) / 0.2e1 + (t67 * mrSges(4,1) - t68 * mrSges(4,2) + Ifges(4,3) * t78 / 0.2e1) * t78 + (t59 * mrSges(5,1) + t57 * mrSges(6,1) - t60 * mrSges(5,2) - t58 * mrSges(6,2) + (Ifges(6,3) / 0.2e1 + Ifges(5,3) / 0.2e1) * t76) * t76 + (t74 * mrSges(5,2) + t65 * mrSges(6,2) - t59 * mrSges(5,3) - t57 * mrSges(6,3) + (Ifges(5,1) / 0.2e1 + Ifges(6,1) / 0.2e1) * t70 + (Ifges(5,5) + Ifges(6,5)) * t76) * t70 + (-t74 * mrSges(5,1) - t65 * mrSges(6,1) + t60 * mrSges(5,3) + t58 * mrSges(6,3) + (Ifges(5,2) / 0.2e1 + Ifges(6,2) / 0.2e1) * t69 + (Ifges(5,6) + Ifges(6,6)) * t76 + (Ifges(5,4) + Ifges(6,4)) * t70) * t69 + ((-t81 * mrSges(3,1) + Ifges(3,2) * t95 / 0.2e1) * t85 + (t81 * mrSges(3,2) + (Ifges(3,4) * t85 + (Ifges(3,1) / 0.2e1 + (qJ(2) * mrSges(4,2) + Ifges(4,1) * t89 / 0.2e1) * t89 + (qJ(2) * mrSges(4,1) - Ifges(4,4) * t89 + Ifges(4,2) * t87 / 0.2e1) * t87) * t84) * qJD(1) + (-t67 * t89 - t68 * t87) * mrSges(4,3) + t78 * (Ifges(4,5) * t89 - Ifges(4,6) * t87)) * t84) * qJD(1) + (Ifges(2,3) / 0.2e1 + (t82 + t83) * qJ(2) * mrSges(3,3)) * t90;
T = t1;
