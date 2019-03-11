% Return the minimum parameter vector for
% S6RPRRRP12
% Use Code from Maple symbolic Code Generation
% 
% Input:
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d1,d3,d4,d5,theta2]';
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
% Datum: 2019-03-09 06:53
% Revision: 8e0af74c1e634ead9bab9e082796ada77f031ee9 (2019-03-08)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function MPV = S6RPRRRP12_convert_par2_MPV_fixb(pkin, m, mrSges, Ifges)

%% Coder Information
%#codegen
%$cgargs {zeros(12,1),zeros(7,1),zeros(7,3),zeros(7,6)}
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RPRRRP12_convert_par2_MPV_fixb: pkin has to be [12x1] (double)');
assert(isreal(m) && all(size(m) == [7 1]), ...
  'S6RPRRRP12_convert_par2_MPV_fixb: m has to be [7x1] (double)'); 
assert(isreal(mrSges) && all(size(mrSges) == [7,3]), ...
  'S6RPRRRP12_convert_par2_MPV_fixb: mrSges has to be [7x3] (double)');
assert(isreal(Ifges) && all(size(Ifges) == [7 6]), ...
  'S6RPRRRP12_convert_par2_MPV_fixb: Ifges has to be [7x6] (double)'); 

%% Symbolic Calculation
% From minimal_parameter_vector_fixb_matlab.m
t199 = (m(5) + m(6));
t190 = m(4) + t199;
t220 = (pkin(9) * t190);
t194 = sin(pkin(7));
t210 = mrSges(4,3) + t220;
t219 = t210 * t194;
t197 = cos(pkin(7));
t218 = t210 * t197;
t205 = (pkin(3) ^ 2);
t185 = (t205 * t199 + Ifges(4,2));
t207 = t185 + (2 * mrSges(4,3) + t220) * pkin(9);
t217 = (-Ifges(6,2) - Ifges(7,3));
t204 = (pkin(4) ^ 2);
t216 = (t204 * m(6) + Ifges(5,2));
t215 = pkin(2) ^ 2 * t190;
t213 = 2 * pkin(11) * mrSges(6,3) - t217;
t212 = 2 * pkin(10) * mrSges(5,3) + t216;
t211 = pkin(11) * m(6) + mrSges(6,3);
t209 = pkin(10) * t199 + mrSges(5,3);
t202 = pkin(10) ^ 2;
t201 = pkin(11) ^ 2;
t198 = cos(pkin(6));
t196 = cos(pkin(12));
t195 = sin(pkin(6));
t193 = sin(pkin(12));
t1 = [Ifges(2,3) + t198 ^ 2 * (t207 * t194 ^ 2 + Ifges(3,3) + t215) + (0.2e1 * (t193 * (-pkin(2) * t218 + Ifges(3,5)) + t196 * (t207 * t197 * t194 + Ifges(3,6))) * t198 + (t196 ^ 2 * (t207 * t197 ^ 2 + Ifges(3,2) + t215) + (0.2e1 * t196 * (pkin(2) * t219 + Ifges(3,4)) + (Ifges(3,1) + t207) * t193) * t193) * t195) * t195; mrSges(2,1); mrSges(2,2); pkin(2) * t190 + mrSges(3,1); mrSges(3,2) - t219; mrSges(3,3) + t218; m(3) + t190; t202 * t199 + Ifges(4,1) - t185 + t212; t209 * pkin(3) + Ifges(4,4); Ifges(4,5); Ifges(4,6); Ifges(4,3) + (t202 + t205) * t199 + t212; pkin(3) * t199 + mrSges(4,1); mrSges(4,2) - t209; t201 * m(6) + Ifges(5,1) + t213 - t216; t211 * pkin(4) + Ifges(5,4); Ifges(5,5); Ifges(5,6); Ifges(5,3) + (t201 + t204) * m(6) + t213; m(6) * pkin(4) + mrSges(5,1); mrSges(5,2) - t211; Ifges(6,1) + Ifges(7,1) + t217; Ifges(6,4) - Ifges(7,5); Ifges(6,5) + Ifges(7,4); Ifges(6,6) - Ifges(7,6); Ifges(6,3) + Ifges(7,2); mrSges(6,1); mrSges(6,2); mrSges(7,1); mrSges(7,2); mrSges(7,3); m(7);];
MPV  = t1;
