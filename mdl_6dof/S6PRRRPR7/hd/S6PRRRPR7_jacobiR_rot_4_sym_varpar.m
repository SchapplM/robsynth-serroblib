% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRRPR7
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [13x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d4,d6,theta1,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:57
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR7_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR7_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRRPR7_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:57:15
% EndTime: 2019-02-22 09:57:15
% DurationCPUTime: 0.15s
% Computational Cost: add. (100->38), mult. (304->88), div. (0->0), fcn. (425->12), ass. (0->45)
t152 = sin(pkin(7));
t153 = sin(pkin(6));
t178 = t152 * t153;
t156 = cos(pkin(6));
t177 = t152 * t156;
t157 = sin(qJ(4));
t176 = t152 * t157;
t160 = cos(qJ(4));
t175 = t152 * t160;
t155 = cos(pkin(7));
t174 = t153 * t155;
t158 = sin(qJ(3));
t173 = t155 * t158;
t161 = cos(qJ(3));
t172 = t155 * t161;
t159 = sin(qJ(2));
t171 = t156 * t159;
t162 = cos(qJ(2));
t170 = t156 * t162;
t169 = t158 * t159;
t168 = t158 * t162;
t167 = t159 * t161;
t166 = t161 * t162;
t165 = t159 * t178;
t151 = sin(pkin(12));
t154 = cos(pkin(12));
t146 = -t151 * t159 + t154 * t170;
t164 = t146 * t155 - t154 * t178;
t148 = -t151 * t170 - t154 * t159;
t163 = t148 * t155 + t151 * t178;
t149 = -t151 * t171 + t154 * t162;
t147 = t151 * t162 + t154 * t171;
t145 = t156 * t155 - t162 * t178;
t144 = (-t155 * t169 + t166) * t153;
t143 = -t148 * t152 + t151 * t174;
t142 = -t146 * t152 - t154 * t174;
t141 = t158 * t177 + (t155 * t168 + t167) * t153;
t140 = t161 * t177 + (t155 * t166 - t169) * t153;
t139 = t148 * t161 - t149 * t173;
t138 = t146 * t161 - t147 * t173;
t137 = t149 * t161 + t163 * t158;
t136 = -t149 * t158 + t163 * t161;
t135 = t147 * t161 + t164 * t158;
t134 = -t147 * t158 + t164 * t161;
t1 = [0, t139 * t160 + t149 * t176, t136 * t160, -t137 * t157 + t143 * t160, 0, 0; 0, t138 * t160 + t147 * t176, t134 * t160, -t135 * t157 + t142 * t160, 0, 0; 0, t144 * t160 + t157 * t165, t140 * t160, -t141 * t157 + t145 * t160, 0, 0; 0, -t139 * t157 + t149 * t175, -t136 * t157, -t137 * t160 - t143 * t157, 0, 0; 0, -t138 * t157 + t147 * t175, -t134 * t157, -t135 * t160 - t142 * t157, 0, 0; 0, -t144 * t157 + t160 * t165, -t140 * t157, -t141 * t160 - t145 * t157, 0, 0; 0, t148 * t158 + t149 * t172, t137, 0, 0, 0; 0, t146 * t158 + t147 * t172, t135, 0, 0, 0; 0 (t155 * t167 + t168) * t153, t141, 0, 0, 0;];
JR_rot  = t1;
