% Return the minimum parameter vector for
% S6RPRRPP4
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [10x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,d1,d3,d4,theta2,theta5]';
% m_mdh [7x1]
%   mass of all robot links (including the base)
% mrSges [7x3]
%  first moment of all robot links (mass times center of mass in body frames)
%  rows: links of the robot (starting with base)
%  columns: x-, y-, z-coordinates
% Ifges [7x6]
%   inertia of all robot links about their respective body frame origins, in body frames
%   rows: links of the robot (starting with base)
%   columns: xx, yy, zz, xy, xz, yz (see inertial_parameters_convert_par1_par2.m)
% 
% Output:
% MPV [27x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 04:41
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPP4_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(10,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [10 1]), ...
  'S6RPRRPP4_convert_par2_MPV_fixb: pkin has to be [10x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPP4_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPP4_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPP4_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t136 = (m(4) + m(5));
t140 = (pkin(3) ^ 2);
t148 = (t140 * m(5) + Ifges(4,2));
t127 = sin(pkin(10));
t129 = cos(pkin(10));
t147 = t127 * t129;
t146 = 2 * pkin(7) * mrSges(4,3) + t148;
t123 = t127 ^ 2;
t124 = t129 ^ 2;
t132 = Ifges(6,2) + Ifges(7,3);
t135 = Ifges(6,1) + Ifges(7,1);
t145 = t123 * t135 + t124 * t132 + Ifges(5,2);
t144 = pkin(8) * m(5) + mrSges(5,3);
t134 = Ifges(6,4) - Ifges(7,5);
t143 = t134 * t147;
t142 = (2 * pkin(8) * mrSges(5,3)) + 0.2e1 * t143 + t145;
t141 = t129 * mrSges(6,1) - t127 * mrSges(6,2);
t139 = pkin(7) ^ 2;
t138 = pkin(8) ^ 2;
t133 = Ifges(6,5) + Ifges(7,4);
t131 = Ifges(6,6) - Ifges(7,6);
t130 = cos(pkin(9));
t128 = sin(pkin(9));
t1 = [Ifges(2,3) + t130 ^ 2 * (Ifges(3,2) + (pkin(2) ^ 2 + t139) * t136 + t146) + (0.2e1 * t130 * Ifges(3,4) + (t139 * t136 + Ifges(3,1) + t146) * t128) * t128; mrSges(2,1); mrSges(2,2); pkin(2) * t136 + mrSges(3,1); mrSges(3,2); pkin(7) * t136 + mrSges(3,3) + mrSges(4,3); m(3) + t136; (t138 * m(5)) + Ifges(4,1) + t142 - t148; t144 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + ((t138 + t140) * m(5)) + t142; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t144; t123 * t132 + t124 * t135 + Ifges(5,1) - 0.4e1 * t143 - t145; Ifges(5,4) + (t124 - t123) * t134 + (-t132 + t135) * t147; -t127 * t131 + t129 * t133 + Ifges(5,5); t127 * t133 + t129 * t131 + Ifges(5,6); 0.2e1 * pkin(4) * t141 + Ifges(7,2) + Ifges(5,3) + Ifges(6,3); mrSges(5,1) + t141; t127 * mrSges(6,1) + t129 * mrSges(6,2) + mrSges(5,2); mrSges(6,3); m(6); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
