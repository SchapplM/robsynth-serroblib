% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6PRRRPR6
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [11x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d2,d3,d4,d6,theta1]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 20:13
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRRPR6_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRRPR6_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6PRRRPR6_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 20:13:35
% EndTime: 2019-02-26 20:13:35
% DurationCPUTime: 0.10s
% Computational Cost: add. (162->46), mult. (461->91), div. (0->0), fcn. (650->12), ass. (0->51)
t152 = sin(pkin(6));
t157 = sin(qJ(3));
t178 = t152 * t157;
t161 = cos(qJ(3));
t177 = t152 * t161;
t162 = cos(qJ(2));
t176 = t152 * t162;
t154 = cos(pkin(6));
t158 = sin(qJ(2));
t175 = t154 * t158;
t174 = t154 * t162;
t156 = sin(qJ(4));
t173 = t156 * t161;
t160 = cos(qJ(4));
t172 = t160 * t161;
t171 = t161 * t162;
t151 = sin(pkin(11));
t153 = cos(pkin(11));
t145 = t151 * t162 + t153 * t175;
t135 = t145 * t161 - t153 * t178;
t144 = t151 * t158 - t153 * t174;
t126 = t135 * t156 - t144 * t160;
t127 = t135 * t160 + t144 * t156;
t155 = sin(qJ(6));
t159 = cos(qJ(6));
t170 = t126 * t159 - t127 * t155;
t169 = t126 * t155 + t127 * t159;
t147 = -t151 * t175 + t153 * t162;
t137 = t147 * t161 + t151 * t178;
t146 = t151 * t174 + t153 * t158;
t128 = t137 * t156 - t146 * t160;
t129 = t137 * t160 + t146 * t156;
t168 = t128 * t159 - t129 * t155;
t167 = t128 * t155 + t129 * t159;
t149 = t154 * t157 + t158 * t177;
t138 = t149 * t156 + t160 * t176;
t139 = t149 * t160 - t156 * t176;
t166 = t138 * t159 - t139 * t155;
t165 = t138 * t155 + t139 * t159;
t164 = -t155 * t160 + t156 * t159;
t163 = t155 * t156 + t159 * t160;
t148 = t154 * t161 - t158 * t178;
t141 = (t156 * t158 + t160 * t171) * t152;
t140 = (t156 * t171 - t158 * t160) * t152;
t136 = -t147 * t157 + t151 * t177;
t134 = -t145 * t157 - t153 * t177;
t133 = -t146 * t172 + t147 * t156;
t132 = -t146 * t173 - t147 * t160;
t131 = -t144 * t172 + t145 * t156;
t130 = -t144 * t173 - t145 * t160;
t1 = [0, t132 * t155 + t133 * t159, t163 * t136, -t168, 0, t168; 0, t130 * t155 + t131 * t159, t163 * t134, -t170, 0, t170; 0, t140 * t155 + t141 * t159, t163 * t148, -t166, 0, t166; 0, t132 * t159 - t133 * t155, t164 * t136, t167, 0, -t167; 0, t130 * t159 - t131 * t155, t164 * t134, t169, 0, -t169; 0, t140 * t159 - t141 * t155, t164 * t148, t165, 0, -t165; 0, t146 * t157, -t137, 0, 0, 0; 0, t144 * t157, -t135, 0, 0, 0; 0, -t157 * t176, -t149, 0, 0, 0;];
JR_rot  = t1;
