% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRPR11
% Use Code from Maple symbolic Code Generation
%
% Rotationsmatrix-Jacobi-Matrix: Differentieller Zusammenhang zwischen
% gestapelter Endeffektor-Rotationsmatrix und verallgemeinerten Koordinaten.
%
%
% Input:
% qJ [6x1]
%   Generalized joint coordinates (joint angles)
% pkin [12x1]
%   kinematic parameters (e.g. lengths of the links)
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d6,theta5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:36
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRPR11_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(12,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRPR11_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [12 1]), ...
  'S6RRRRPR11_jacobiR_rot_6_sym_varpar: pkin has to be [12x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:36:21
% EndTime: 2019-02-26 22:36:21
% DurationCPUTime: 0.11s
% Computational Cost: add. (194->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t157 = cos(pkin(6));
t159 = sin(qJ(2));
t163 = cos(qJ(1));
t165 = t163 * t159;
t160 = sin(qJ(1));
t162 = cos(qJ(2));
t167 = t160 * t162;
t148 = t157 * t165 + t167;
t158 = sin(qJ(3));
t161 = cos(qJ(3));
t156 = sin(pkin(6));
t169 = t156 * t163;
t142 = -t148 * t161 + t158 * t169;
t164 = t163 * t162;
t168 = t160 * t159;
t147 = -t157 * t164 + t168;
t155 = qJ(4) + pkin(12) + qJ(6);
t153 = sin(t155);
t154 = cos(t155);
t134 = t142 * t153 + t147 * t154;
t135 = t142 * t154 - t147 * t153;
t174 = t153 * t161;
t173 = t154 * t161;
t172 = t156 * t158;
t171 = t156 * t161;
t170 = t156 * t162;
t166 = t161 * t162;
t140 = -t148 * t158 - t161 * t169;
t150 = -t157 * t168 + t164;
t149 = t157 * t167 + t165;
t146 = t157 * t158 + t159 * t171;
t145 = t157 * t161 - t159 * t172;
t144 = t150 * t161 + t160 * t172;
t143 = t150 * t158 - t160 * t171;
t139 = -t146 * t154 + t153 * t170;
t138 = -t146 * t153 - t154 * t170;
t137 = t144 * t154 + t149 * t153;
t136 = -t144 * t153 + t149 * t154;
t1 = [t135, -t149 * t173 + t150 * t153, -t143 * t154, t136, 0, t136; t137, -t147 * t173 + t148 * t153, t140 * t154, t134, 0, t134; 0 (t153 * t159 + t154 * t166) * t156, t145 * t154, t138, 0, t138; -t134, t149 * t174 + t150 * t154, t143 * t153, -t137, 0, -t137; t136, t147 * t174 + t148 * t154, -t140 * t153, t135, 0, t135; 0 (-t153 * t166 + t154 * t159) * t156, -t145 * t153, t139, 0, t139; t140, -t149 * t158, t144, 0, 0, 0; t143, -t147 * t158, -t142, 0, 0, 0; 0, t158 * t170, t146, 0, 0, 0;];
JR_rot  = t1;
