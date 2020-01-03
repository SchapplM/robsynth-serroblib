% Return the minimum parameter vector for
% S5RPRRR14
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-12-31 19:22
% Revision: 9bd3e9fa678258af3b32f1bcc8622e39ff85504d (2019-12-30)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S5RPRRR14_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(11,1),zeros(6,1),zeros(6,3),zeros(6,6)}
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S5RPRRR14_convert_par2_MPV_fixb: pkin has to be [11x1] (double)');
assert(isreal(m) && all(size(m) == [6 1]), ...
  'S5RPRRR14_convert_par2_MPV_fixb: m has to be [6x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [6,3]), ...
  'S5RPRRR14_convert_par2_MPV_fixb: mrSges has to be [6x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [6 6]), ...
  'S5RPRRR14_convert_par2_MPV_fixb: Ifges has to be [6x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t170 = (m(5) + m(6));
t161 = m(4) + t170;
t190 = (pkin(8) * t161);
t165 = sin(pkin(6));
t181 = mrSges(4,3) + t190;
t189 = t181 * t165;
t168 = cos(pkin(6));
t188 = t181 * t168;
t176 = (pkin(3) ^ 2);
t156 = (t176 * t170 + Ifges(4,2));
t178 = t156 + (2 * mrSges(4,3) + t190) * pkin(8);
t175 = (pkin(4) ^ 2);
t187 = (t175 * m(6) + Ifges(5,2));
t186 = 2 * pkin(10) * mrSges(6,3) + Ifges(6,2);
t185 = pkin(2) ^ 2 * t161;
t183 = 2 * pkin(9) * mrSges(5,3) + t187;
t182 = pkin(10) * m(6) + mrSges(6,3);
t180 = pkin(9) * t170 + mrSges(5,3);
t173 = pkin(9) ^ 2;
t172 = pkin(10) ^ 2;
t169 = cos(pkin(5));
t167 = cos(pkin(11));
t166 = sin(pkin(5));
t164 = sin(pkin(11));
t1 = [Ifges(2,3) + t169 ^ 2 * (t178 * t165 ^ 2 + Ifges(3,3) + t185) + (0.2e1 * (t164 * (-pkin(2) * t188 + Ifges(3,5)) + t167 * (t178 * t168 * t165 + Ifges(3,6))) * t169 + (t167 ^ 2 * (t178 * t168 ^ 2 + Ifges(3,2) + t185) + (0.2e1 * t167 * (pkin(2) * t189 + Ifges(3,4)) + (Ifges(3,1) + t178) * t164) * t164) * t166) * t166; mrSges(2,1); mrSges(2,2); pkin(2) * t161 + mrSges(3,1); mrSges(3,2) - t189; mrSges(3,3) + t188; m(3) + t161; t173 * t170 + Ifges(4,1) - t156 + t183; t180 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t173 + t176) * t170 + t183; pkin(3) * t170 + mrSges(4,1); mrSges(4,2) - t180; m(6) * t172 + Ifges(5,1) + t186 - t187; t182 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t172 + t175) * m(6) + t186; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t182; Ifges(6,1) - Ifges(6,2); Ifges(6,4); Ifges(6,5); Ifges(6,6); Ifges(6,3); mrSges(6,1); mrSges(6,2);];
MPV = t1;
