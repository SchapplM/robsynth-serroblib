% Return the minimum parameter vector for
% S6RPRRPR11
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d6,theta2,theta5]';
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
% MPV [32x1]
%   base parameter vector (minimal parameter vector)

% Quelle: HybrDyn-Toolbox
% Datum: 2019-03-09 05:47
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRPR11_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(13,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6RPRRPR11_convert_par2_MPV_fixb: pkin has to be [13x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRPR11_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRPR11_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRPR11_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t201 = m(4) + m(5);
t224 = (pkin(9) * t201);
t195 = sin(pkin(7));
t213 = mrSges(4,3) + t224;
t223 = t213 * t195;
t199 = cos(pkin(7));
t222 = t213 * t199;
t208 = (pkin(3) ^ 2);
t185 = (t208 * m(5) + Ifges(4,2));
t210 = t185 + (2 * mrSges(4,3) + t224) * pkin(9);
t221 = 2 * pkin(11) * mrSges(7,3) + Ifges(7,2);
t193 = sin(pkin(13));
t197 = cos(pkin(13));
t220 = t193 * t197;
t219 = pkin(2) ^ 2 * t201;
t207 = (pkin(5) ^ 2);
t217 = (t207 * m(7) + Ifges(5,2) + Ifges(6,3));
t216 = Ifges(6,4) * t220;
t215 = pkin(10) * m(5) + mrSges(5,3);
t214 = -pkin(11) * m(7) - mrSges(7,3);
t212 = 2 * pkin(10) * mrSges(5,3) + t217;
t205 = pkin(10) ^ 2;
t204 = pkin(11) ^ 2;
t200 = cos(pkin(6));
t198 = cos(pkin(12));
t196 = sin(pkin(6));
t194 = sin(pkin(12));
t189 = t197 ^ 2;
t186 = t193 ^ 2;
t184 = t214 * pkin(5) + Ifges(6,5);
t183 = m(7) * t204 + Ifges(6,1) + t221;
t182 = Ifges(6,2) + (t204 + t207) * m(7) + t221;
t1 = [Ifges(2,3) + t200 ^ 2 * (t210 * t195 ^ 2 + Ifges(3,3) + t219) + (0.2e1 * (t194 * (-pkin(2) * t222 + Ifges(3,5)) + t198 * (t210 * t199 * t195 + Ifges(3,6))) * t200 + (t198 ^ 2 * (t210 * t199 ^ 2 + Ifges(3,2) + t219) + (0.2e1 * t198 * (pkin(2) * t223 + Ifges(3,4)) + (Ifges(3,1) + t210) * t194) * t194) * t196) * t196; mrSges(2,1); mrSges(2,2); pkin(2) * t201 + mrSges(3,1); mrSges(3,2) - t223; mrSges(3,3) + t222; m(3) + t201; t205 * m(5) + Ifges(4,1) - t185 + t212; t215 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t205 + t208) * m(5) + t212; m(5) * pkin(3) + mrSges(4,1); mrSges(4,2) - t215; t186 * t182 + t189 * t183 + Ifges(5,1) - 0.2e1 * t216 - t217; t193 * Ifges(6,6) - t197 * t184 + Ifges(5,4); Ifges(5,5) + (t189 - t186) * Ifges(6,4) + (-t182 + t183) * t220; -t197 * Ifges(6,6) - t193 * t184 + Ifges(5,6); t189 * t182 + t186 * t183 + Ifges(5,3) + 0.2e1 * t216; mrSges(5,1); mrSges(5,2); m(7) * pkin(5) + mrSges(6,1); mrSges(6,2); mrSges(6,3) - t214; m(6) + m(7); Ifges(7,1) - Ifges(7,2); Ifges(7,4); Ifges(7,5); Ifges(7,6); Ifges(7,3); mrSges(7,1); mrSges(7,2);];
MPV  = t1;
