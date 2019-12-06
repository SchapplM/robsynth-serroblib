% Return the minimum parameter vector for
% S5RPRRR2
% Use Code from Maple symbolic Code Generation
% 
% Input:
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
% MPV [28x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-12-05 18:13
% Revision: 77da58f92bca3eff71542919beafa37024070d86 (2019-12-05)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRRR2_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(9,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [9 1]), ...
  'S5RPRRR2_convert_par2_MPV_fixb: pkin has to be [9x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR2_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR2_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR2_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t106 = -pkin(8) * m(6) - mrSges(6,3);
t86 = (mrSges(5,3) - t106);
t94 = (m(5) + m(6));
t105 = -pkin(7) * t94 - t86;
t100 = (pkin(3) ^ 2);
t97 = (pkin(7) ^ 2);
t104 = (Ifges(4,2) + (t97 + t100) * t94);
t96 = (pkin(8) ^ 2);
t99 = (pkin(4) ^ 2);
t103 = (Ifges(5,2) + (t96 + t99) * m(6));
t90 = (m(4) + t94);
t102 = (mrSges(4,3) - t105);
t101 = 2 * pkin(8) * mrSges(6,3) + 2 * pkin(6) * t102 + 2 * pkin(7) * t86 + Ifges(6,2) + t103 + t104;
t98 = pkin(6) ^ 2;
t93 = cos(pkin(9));
t92 = sin(pkin(9));
t1 = [Ifges(2,3) + t93 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t98) * t90 + t101) + (0.2e1 * t93 * Ifges(3,4) + (t98 * t90 + Ifges(3,1) + t101) * t92) * t92; mrSges(2,1); mrSges(2,2); pkin(2) * t90 + mrSges(3,1); mrSges(3,2); pkin(6) * t90 + mrSges(3,3) + t102; m(3) + t90; t97 * t94 + Ifges(4,1) - t104; Ifges(4,4); t105 * pkin(3) + Ifges(4,5); Ifges(4,6); t100 * t94 + Ifges(4,3); pkin(3) * t94 + mrSges(4,1); mrSges(4,2); m(6) * t96 + Ifges(5,1) - t103; Ifges(5,4); t106 * pkin(4) + Ifges(5,5); Ifges(5,6); t99 * m(6) + Ifges(5,3); m(6) * pkin(4) + mrSges(5,1); mrSges(5,2); Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
