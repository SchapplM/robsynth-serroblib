% Rotatorische Teilmatrix der Rotationsmatrix-Jacobi-Matrix für Segment Nr. 5 (0=Basis) von
% S6RRRPPR9
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
%   pkin=[a2,a3,a4,a5,a6,alpha2,d1,d2,d3,d6,theta4]';
%
% Output:
% JR_rot [9x6]
%   Jacobi-Matrix der Endeffektor-Rotationsmatrix

% Quelle: HybrDyn-Toolbox
% Datum: 2019-02-26 22:08
% Revision: d75aae1ac561373cd3be920984c3df29a1c2ecc8 (2019-02-26)
% Moritz Schappler, moritz.schappler@imes.uni-hannover.de
% (C) Institut für Mechatronische Systeme, Universität Hannover

function JR_rot = S6RRRPPR9_jacobiR_rot_5_sym_varpar(qJ, ...
  pkin)
%% Coder Information
%#codegen
%$cgargs {zeros(6,1),zeros(11,1)}
assert(isreal(qJ) && all(size(qJ) == [6 1]), ...
  'S6RRRPPR9_jacobiR_rot_5_sym_varpar: qJ has to be [6x1] (double)');
assert(isreal(pkin) && all(size(pkin) == [11 1]), ...
  'S6RRRPPR9_jacobiR_rot_5_sym_varpar: pkin has to be [11x1] (double)');

%% Symbolic Calculation
% From jacobiR_rot_5_floatb_twist_matlab.m
% OptimizationMode: 2
% StartTime: 2019-02-26 22:08:10
% EndTime: 2019-02-26 22:08:10
% DurationCPUTime: 0.07s
% Computational Cost: add. (51->24), mult. (163->57), div. (0->0), fcn. (238->10), ass. (0->30)
t136 = sin(pkin(11));
t143 = cos(qJ(3));
t155 = t136 * t143;
t137 = sin(pkin(6));
t140 = sin(qJ(3));
t154 = t137 * t140;
t153 = t137 * t143;
t145 = cos(qJ(1));
t152 = t137 * t145;
t138 = cos(pkin(11));
t151 = t138 * t143;
t141 = sin(qJ(2));
t142 = sin(qJ(1));
t150 = t142 * t141;
t144 = cos(qJ(2));
t149 = t142 * t144;
t148 = t143 * t144;
t147 = t145 * t141;
t146 = t145 * t144;
t139 = cos(pkin(6));
t132 = t139 * t147 + t149;
t126 = -t132 * t140 - t143 * t152;
t127 = -t132 * t143 + t140 * t152;
t134 = -t139 * t150 + t146;
t133 = t139 * t149 + t147;
t131 = t139 * t146 - t150;
t130 = t139 * t143 - t141 * t154;
t129 = t134 * t143 + t142 * t154;
t128 = t134 * t140 - t142 * t153;
t1 = [t127 * t138 + t131 * t136, -t133 * t151 + t134 * t136, -t128 * t138, 0, 0, 0; t129 * t138 + t133 * t136, t131 * t151 + t132 * t136, t126 * t138, 0, 0, 0; 0 (t136 * t141 + t138 * t148) * t137, t130 * t138, 0, 0, 0; t126, -t133 * t140, t129, 0, 0, 0; t128, t131 * t140, -t127, 0, 0, 0; 0, t144 * t154, t139 * t140 + t141 * t153, 0, 0, 0; t127 * t136 - t131 * t138, -t133 * t155 - t134 * t138, -t128 * t136, 0, 0, 0; t129 * t136 - t133 * t138, t131 * t155 - t132 * t138, t126 * t136, 0, 0, 0; 0 (t136 * t148 - t138 * t141) * t137, t130 * t136, 0, 0, 0;];
JR_rot  = t1;
