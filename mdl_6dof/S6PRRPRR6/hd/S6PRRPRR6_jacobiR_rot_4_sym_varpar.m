% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 4 (0=Basis) von
% S6PRRPRR6
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,alpha3,d2,d3,d5,d6,theta1,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-22 09:51
% Revision: 2b76964ad985d937eecd005a1a368749e6b3dc4d (2019-02-18)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6PRRPRR6_jacobiR_rot_4_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(13,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6PRRPRR6_jacobiR_rot_4_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [13 1]), ...
  'S6PRRPRR6_jacobiR_rot_4_sym_varpar: pkin has to be [13x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_4_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-22 09:50:54
% EndTime: 2019-02-22 09:50:54
% DurationCPUTime: 0.13s
% Computational Cost: add. (64->29), mult. (200->69), div. (0->0), fcn. (279->12), ass. (0->38)
t139 = sin(pkin(13));
t141 = sin(pkin(7));
t165 = t139 * t141;
t142 = sin(pkin(6));
t164 = t141 * t142;
t143 = cos(pkin(13));
t163 = t141 * t143;
t146 = cos(pkin(6));
t162 = t141 * t146;
t145 = cos(pkin(7));
t147 = sin(qJ(3));
t161 = t145 * t147;
t149 = cos(qJ(3));
t160 = t145 * t149;
t148 = sin(qJ(2));
t159 = t146 * t148;
t150 = cos(qJ(2));
t158 = t146 * t150;
t157 = t147 * t148;
t156 = t147 * t150;
t155 = t148 * t149;
t154 = t149 * t150;
t153 = t148 * t164;
t140 = sin(pkin(12));
t144 = cos(pkin(12));
t134 = -t140 * t148 + t144 * t158;
t152 = t134 * t145 - t144 * t164;
t136 = -t140 * t158 - t144 * t148;
t151 = t136 * t145 + t140 * t164;
t137 = -t140 * t159 + t144 * t150;
t135 = t140 * t150 + t144 * t159;
t133 = (-t145 * t157 + t154) * t142;
t132 = t149 * t162 + (t145 * t154 - t157) * t142;
t131 = t136 * t149 - t137 * t161;
t130 = t134 * t149 - t135 * t161;
t129 = -t137 * t147 + t151 * t149;
t128 = -t135 * t147 + t152 * t149;
t1 = [0, t131 * t143 + t137 * t165, t129 * t143, 0, 0, 0; 0, t130 * t143 + t135 * t165, t128 * t143, 0, 0, 0; 0, t133 * t143 + t139 * t153, t132 * t143, 0, 0, 0; 0, -t131 * t139 + t137 * t163, -t129 * t139, 0, 0, 0; 0, -t130 * t139 + t135 * t163, -t128 * t139, 0, 0, 0; 0, -t133 * t139 + t143 * t153, -t132 * t139, 0, 0, 0; 0, t136 * t147 + t137 * t160, t137 * t149 + t151 * t147, 0, 0, 0; 0, t134 * t147 + t135 * t160, t135 * t149 + t152 * t147, 0, 0, 0; 0 (t145 * t155 + t156) * t142, t147 * t162 + (t145 * t156 + t155) * t142, 0, 0, 0;];
JR_rot  = t1;
