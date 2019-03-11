% Return the minimum parameter vector for
% S6PRRPRR3
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
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
% Datum: 2019-03-08 22:09
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6PRRPRR3_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR3_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6PRRPRR3_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6PRRPRR3_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6PRRPRR3_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t185 = (m(6) + m(7));
t188 = (pkin(10) ^ 2);
t189 = (pkin(5) ^ 2);
t199 = (t189 * m(7) + Ifges(6,2));
t195 = 2 * pkin(10) * mrSges(6,3) + t199;
t172 = t188 * t185 + Ifges(5,1) + t195;
t190 = (pkin(4) ^ 2);
t176 = t190 * t185 + Ifges(5,2);
t201 = t172 - t176;
t200 = (m(4) * pkin(9));
t198 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t182 = sin(pkin(13));
t184 = cos(pkin(13));
t197 = t182 * t184;
t177 = t182 ^ 2;
t179 = t184 ^ 2;
t196 = t179 - t177;
t194 = pkin(11) * m(7) + mrSges(7,3);
t192 = pkin(10) * t185 + mrSges(6,3);
t173 = t192 * pkin(4) + Ifges(5,4);
t193 = t173 * t197;
t174 = mrSges(5,2) - t192;
t175 = pkin(4) * t185 + mrSges(5,1);
t191 = -t182 * t174 + t184 * t175;
t187 = pkin(11) ^ 2;
t183 = sin(pkin(7));
t1 = [m(2) + m(3) + m(4); pkin(2) ^ 2 * m(4) + Ifges(3,3) + (0.2e1 * t193 + t177 * t172 + t179 * t176 + Ifges(4,2) + ((2 * mrSges(4,3) + t200) * pkin(9))) * t183 ^ 2; m(4) * pkin(2) + mrSges(3,1); mrSges(3,2) + (-mrSges(4,3) - t200) * t183; t201 * t196 + Ifges(4,1) - Ifges(4,2) - 0.4e1 * t193; t196 * t173 + t201 * t197 + Ifges(4,4); t184 * Ifges(5,5) - t182 * Ifges(5,6) + Ifges(4,5); t182 * Ifges(5,5) + t184 * Ifges(5,6) + Ifges(4,6); Ifges(4,3) + Ifges(5,3) + ((t188 + t190) * t185) + 0.2e1 * t191 * pkin(3) + t195; mrSges(4,1) + t191; t184 * t174 + t182 * t175 + mrSges(4,2); mrSges(5,3); m(5) + t185; m(7) * t187 + Ifges(6,1) + t198 - t199; t194 * pkin(5) + Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3) + (t187 + t189) * m(7) + t198; m(7) * pkin(5) + mrSges(6,1); mrSges(6,2) - t194; Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
