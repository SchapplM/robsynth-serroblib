% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 6 (0=Basis) von
% S6RRRRRP9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d4,d5]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 12:31
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRRRP9_jacobiR_rot_6_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRRRP9_jacobiR_rot_6_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRRRP9_jacobiR_rot_6_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_6_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 12:30:54
% EndTime: 2019-02-22 12:30:54
% DurationCPUTime: 0.11s
% Computational Cost: add. (144->32), mult. (275->62), div. (0->0), fcn. (402->10), ass. (0->39)
t158 = cos(pkin(6));
t160 = sin(qJ(2));
t164 = cos(qJ(1));
t166 = t164 * t160;
t161 = sin(qJ(1));
t163 = cos(qJ(2));
t168 = t161 * t163;
t149 = t158 * t166 + t168;
t159 = sin(qJ(3));
t162 = cos(qJ(3));
t157 = sin(pkin(6));
t170 = t157 * t164;
t143 = -t149 * t162 + t159 * t170;
t165 = t164 * t163;
t169 = t161 * t160;
t148 = -t158 * t165 + t169;
t156 = qJ(4) + qJ(5);
t154 = sin(t156);
t155 = cos(t156);
t135 = t143 * t154 + t148 * t155;
t136 = t143 * t155 - t148 * t154;
t175 = t154 * t162;
t174 = t155 * t162;
t173 = t157 * t159;
t172 = t157 * t162;
t171 = t157 * t163;
t167 = t162 * t163;
t141 = -t149 * t159 - t162 * t170;
t151 = -t158 * t169 + t165;
t150 = t158 * t168 + t166;
t147 = t158 * t159 + t160 * t172;
t146 = t158 * t162 - t160 * t173;
t145 = t151 * t162 + t161 * t173;
t144 = t151 * t159 - t161 * t172;
t140 = -t147 * t155 + t154 * t171;
t139 = -t147 * t154 - t155 * t171;
t138 = t145 * t155 + t150 * t154;
t137 = -t145 * t154 + t150 * t155;
t1 = [t136, -t150 * t174 + t151 * t154, -t144 * t155, t137, t137, 0; t138, -t148 * t174 + t149 * t154, t141 * t155, t135, t135, 0; 0 (t154 * t160 + t155 * t167) * t157, t146 * t155, t139, t139, 0; -t135, t150 * t175 + t151 * t155, t144 * t154, -t138, -t138, 0; t137, t148 * t175 + t149 * t155, -t141 * t154, t136, t136, 0; 0 (-t154 * t167 + t155 * t160) * t157, -t146 * t154, t140, t140, 0; t141, -t150 * t159, t145, 0, 0, 0; t144, -t148 * t159, -t143, 0, 0, 0; 0, t159 * t171, t147, 0, 0, 0;];
JR_rot  = t1;
